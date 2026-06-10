from pathlib import Path

import numpy as np
import rasterio as ras
import xarray as xr
from dask.distributed import Client
from rasterio.transform import from_bounds

from geodata import XarrayMask
from geodata.datasets import load_dataset
from geodata.model.wind import WindInterpolationModel


def _create_saved_mask_from_output_grid(
    output: xr.DataArray,
    mask_dir: Path,
    name: str = "wind_xmask",
) -> None:
    x = output["x"].values
    y = output["y"].values

    lon = np.sort(np.asarray(x, dtype=float))
    lat = np.sort(np.asarray(y, dtype=float))
    lon_step = float(np.abs(lon[1] - lon[0]))
    lat_step = float(np.abs(lat[1] - lat[0]))
    west = float(lon.min() - lon_step / 2)
    east = float(lon.max() + lon_step / 2)
    south = float(lat.min() - lat_step / 2)
    north = float(lat.max() + lat_step / 2)

    nlon_hi = len(lon) * 2
    nlat_hi = len(lat) * 2
    transform = from_bounds(west, south, east, north, nlon_hi, nlat_hi)

    arr = np.zeros((nlat_hi, nlon_hi), dtype=np.uint8)
    arr[nlat_hi // 4 : 3 * nlat_hi // 4, nlon_hi // 4 : 3 * nlon_hi // 4] = 1

    source_tif = mask_dir / "source.tif"
    with ras.open(
        str(source_tif),
        "w",
        driver="GTiff",
        height=arr.shape[0],
        width=arr.shape[1],
        count=1,
        dtype=arr.dtype,
        compress="lzw",
        crs="+proj=latlong",
        transform=transform,
    ) as dst:
        dst.write(arr, 1)

    from geodata import Mask

    mask = Mask(name=name, mask_dir=str(mask_dir))
    mask.add_layer(str(source_tif), layer_name="source")
    mask.merge_layer(show_raster=False)
    mask.save_mask()


def test_wind_estimate_with_xarray_mask_offline(tmp_path):
    years = slice(2016, 2016)
    months = slice(1, 1)

    with Client(processes=True, threads_per_worker=1):
        ds_cls = load_dataset("wind_3d_hourly_test")
        ds = ds_cls(years=years, months=months)
        assert ds.downloaded, "Wind fixture NetCDF should be present"

        model = WindInterpolationModel(ds)
        model.prepare(force=True)

        base = model.estimate(years=years, months=months, height=12)
        assert isinstance(base, xr.DataArray)

        mask_dir = tmp_path / "masks"
        mask_dir.mkdir(parents=True, exist_ok=True)
        mask_name = "wind_xmask"
        _create_saved_mask_from_output_grid(base, mask_dir, name=mask_name)

        base_ds = base.to_dataset(name=base.name or "value")
        xmask = XarrayMask.from_name(mask_name, grid=base_ds, mask_dir=str(mask_dir))
        masked = xmask.apply(
            base_ds,
            mode="where",
            include_area=True,
        )

        assert isinstance(masked, dict)
        assert set(masked.keys()) == {"merged_mask"}

        merged = masked["merged_mask"]
        assert "area" in merged
        value_vars = [v for v in merged.data_vars if v not in {"area"}]
        assert len(value_vars) == 1
        var = value_vars[0]

        attached = xmask.attach(base, include_area=False)["merged_mask"]
        valid = attached["mask"] > 0
        expected = attached[var].where(valid)
        xr.testing.assert_allclose(merged[var], expected)
