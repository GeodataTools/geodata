from pathlib import Path
from typing import Any, cast

import numpy as np
import shapely.geometry
import xarray as xr
import rasterio as ras
from rasterio.transform import from_bounds

from geodata import Mask, XarrayMask
from geodata.cutout import Cutout, ds_reformat_index


def _build_minimal_cutout() -> Cutout:
    cutout = Cutout.__new__(Cutout)
    cutout.name = "xarray-mask-test"
    cutout.meta = xr.Dataset(
        coords={
            "x": np.array([100.0, 100.25, 100.5, 100.75]),
            "y": np.array([30.75, 30.5, 30.25, 30.0]),
            "year": [2016],
            "month": [1],
        }
    )
    cutout.merged_mask = None
    cutout.shape_mask = None
    cutout.area = None
    cutout.prepared = True
    cutout.empty = False
    cutout.cutout_dir = Path(".")
    return cutout


def _create_saved_mask(cutout: Cutout, mask_dir: Path, name: str = "xarray_test_mask") -> None:
    assert cutout.meta is not None
    xr_meta = ds_reformat_index(cast(Any, cutout.meta))
    lon = xr_meta["lon"].values
    lat = xr_meta["lat"].values

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
    layer_path = mask_dir / "source.tif"
    with ras.open(
        str(layer_path),
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

    mask = Mask(name=name, mask_dir=str(mask_dir))
    mask.add_layer(str(layer_path), layer_name="source")
    mask.merge_layer(show_raster=False)
    shape = shapely.geometry.box(west, south, (west + east) / 2, (south + north) / 2)
    mask.extract_shapes({"region_a": shape}, show_raster=False)
    mask.save_mask()


def _sample_dataset_from_cutout(cutout: Cutout) -> xr.Dataset:
    assert cutout.meta is not None
    y = cutout.meta["y"].values
    x = cutout.meta["x"].values
    t = np.array(["2016-01-01T00:00:00", "2016-01-01T01:00:00"], dtype="datetime64[ns]")
    vals = np.arange(len(t) * len(y) * len(x), dtype=np.float32).reshape(
        len(t), len(y), len(x)
    )
    return xr.Dataset({"signal": (("time", "y", "x"), vals)}, coords={"time": t, "y": y, "x": x})


def test_xarraymask_attach_matches_legacy_contract(tmp_path, monkeypatch):
    cutout = _build_minimal_cutout()
    mask_dir = tmp_path / "masks"
    mask_dir.mkdir(parents=True, exist_ok=True)
    mask_name = "xarray_test_mask"
    _create_saved_mask(cutout, mask_dir, name=mask_name)

    monkeypatch.setattr("geodata.cutout.config.MASK_DIR", str(mask_dir))
    cutout.add_mask(mask_name)
    cutout.add_grid_area()

    ds = _sample_dataset_from_cutout(cutout)
    legacy = cutout.mask(ds)

    assert cutout.meta is not None
    xmask = XarrayMask.from_name(mask_name, grid=cutout.meta, mask_dir=str(mask_dir))
    attached = xmask.attach(ds, include_area=True)

    assert set(attached.keys()) == set(legacy.keys())
    for key in attached:
        xr.testing.assert_allclose(attached[key]["mask"], legacy[key]["mask"])
        xr.testing.assert_allclose(attached[key]["area"], legacy[key]["area"])
        xr.testing.assert_allclose(attached[key]["signal"], legacy[key]["signal"])


def test_xarraymask_apply_where_and_multiply(tmp_path):
    cutout = _build_minimal_cutout()
    mask_dir = tmp_path / "masks"
    mask_dir.mkdir(parents=True, exist_ok=True)
    mask_name = "xarray_test_mask"
    _create_saved_mask(cutout, mask_dir, name=mask_name)

    ds = _sample_dataset_from_cutout(cutout)
    assert cutout.meta is not None
    xmask = XarrayMask.from_name(mask_name, grid=cutout.meta, mask_dir=str(mask_dir))

    attached = xmask.attach(ds, include_area=False)
    merged_mask = attached["merged_mask"]["mask"]

    where_out = xmask.apply(ds, mode="where", include_area=True)["merged_mask"]
    multiply_out = xmask.apply(ds, mode="multiply", include_area=False)["merged_mask"]

    valid = merged_mask > 0
    expected_where = attached["merged_mask"]["signal"].where(valid)
    expected_multiply = attached["merged_mask"]["signal"] * valid

    xr.testing.assert_allclose(where_out["signal"], expected_where)
    xr.testing.assert_allclose(multiply_out["signal"], expected_multiply)
    assert "area" in where_out
    assert "area" not in multiply_out
