import uuid
from pathlib import Path
from typing import Any, cast

import numpy as np
import shapely.geometry
import xarray as xr
from rasterio.transform import from_bounds

from geodata.cutout import Cutout, calc_grid_area, coarsen, ds_reformat_index
from geodata.datasets import load_dataset
from geodata.mask import Mask, save_raster


def _build_cutout(tmp_path: Path) -> Cutout:
    dataset_cls = load_dataset("wind_solar_hourly_test")
    dataset = dataset_cls(years=slice(2016, 2016), months=slice(1, 1), testing=True)
    assert dataset.downloaded, "Fixture NetCDF should be present"

    with xr.open_dataset(dataset.catalog[0].path, engine="h5netcdf") as opened:
        if "x" in opened.coords and "y" in opened.coords:
            xvals = opened["x"].values
            yvals = opened["y"].values
        else:
            xvals = opened["longitude"].values
            yvals = opened["latitude"].values

    # Use a lightweight Cutout instance that still exercises legacy methods
    # (add_mask, add_grid_area, mask) without invoking dataset preparation.
    cutout = Cutout.__new__(Cutout)
    cutout.name = f"legacy-mask-test-{uuid.uuid4().hex[:8]}"
    cutout.meta = xr.Dataset(
        coords={
            "x": xvals,
            "y": yvals,
            "year": [2016],
            "month": [1],
        }
    )
    cutout.merged_mask = None
    cutout.shape_mask = None
    cutout.area = None
    cutout.prepared = True
    cutout.empty = False
    cutout.cutout_dir = tmp_path / "cutouts"
    return cutout


def _create_and_save_mask(cutout: Cutout, mask_dir: Path, name: str = "legacy_test_mask") -> None:
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

    raster = np.zeros((nlat_hi, nlon_hi), dtype=np.uint8)
    # Non-trivial pattern so coarsening does real work.
    raster[nlat_hi // 4 : 3 * nlat_hi // 4, nlon_hi // 6 : 5 * nlon_hi // 6] = 1

    layer_path = mask_dir / "source_layer.tif"
    save_raster(raster, transform, str(layer_path))

    mask = Mask(name=name, mask_dir=str(mask_dir))
    mask.add_layer(str(layer_path), layer_name="source")
    mask.merge_layer(show_raster=False)

    centroid_lon = float(np.mean([west, east]))
    centroid_lat = float(np.mean([south, north]))
    shape = shapely.geometry.box(
        west,
        south,
        centroid_lon,
        centroid_lat,
    )
    mask.extract_shapes({"region_a": shape}, show_raster=False)
    mask.save_mask()


def test_legacy_mask_workflow_contract_offline(tmp_path, monkeypatch):
    cutout = _build_cutout(tmp_path)
    mask_dir = tmp_path / "masks"
    mask_dir.mkdir(parents=True, exist_ok=True)
    mask_name = "legacy_test_mask"
    _create_and_save_mask(cutout, mask_dir, name=mask_name)

    monkeypatch.setattr("geodata.cutout.config.MASK_DIR", str(mask_dir))
    cutout.add_mask(mask_name)
    cutout.add_grid_area()

    time = np.array(["2016-01-01T00:00:00", "2016-01-01T01:00:00"], dtype="datetime64[ns]")
    y = cutout.coords["y"].values
    x = cutout.coords["x"].values
    payload = np.arange(len(time) * len(y) * len(x), dtype=np.float32).reshape(
        len(time), len(y), len(x)
    )
    ds = xr.Dataset(
        {"signal": (("time", "y", "x"), payload)},
        coords={"time": time, "y": y, "x": x},
    )

    masked = cutout.mask(ds)

    assert set(masked.keys()) == {"merged_mask", "region_a"}
    merged = masked["merged_mask"]
    assert isinstance(merged, xr.Dataset)
    assert {"signal", "mask", "area"}.issubset(set(merged.data_vars))
    assert tuple(merged["signal"].dims) == ("time", "lat", "lon")
    assert tuple(merged["mask"].dims) == ("lat", "lon")
    assert tuple(merged["area"].dims) == ("lat", "lon")


def test_legacy_add_mask_coarsen_parity_offline(tmp_path, monkeypatch):
    cutout = _build_cutout(tmp_path)
    mask_dir = tmp_path / "masks"
    mask_dir.mkdir(parents=True, exist_ok=True)
    mask_name = "legacy_test_mask"
    _create_and_save_mask(cutout, mask_dir, name=mask_name)

    monkeypatch.setattr("geodata.cutout.config.MASK_DIR", str(mask_dir))
    cutout.add_mask(mask_name, shape_mask=False)

    mask = Mask.from_name(mask_name, mask_dir=str(mask_dir))
    assert cutout.meta is not None
    expected = coarsen(
        cast(Any, mask.load_merged_xr()),
        cast(Any, ds_reformat_index(cast(Any, cutout.meta))),
    )

    assert cutout.merged_mask is not None
    np.testing.assert_allclose(cutout.merged_mask.values, expected.values)
    assert cutout.merged_mask.shape == expected.shape


def test_legacy_add_grid_area_sanity_offline(tmp_path):
    cutout = _build_cutout(tmp_path)
    cutout.add_grid_area()

    assert cutout.area is not None
    area = cutout.area["area"].values
    assert np.all(np.isfinite(area))
    assert np.all(area > 0)

    # Area should be constant across longitude for a given latitude row.
    row_std = area.std(axis=1)
    assert np.allclose(row_std, 0.0, atol=1e-6)

    assert cutout.meta is not None
    xr_ds = ds_reformat_index(cast(Any, cutout.meta))
    lat = xr_ds.lat.values
    lon = xr_ds.lon.values
    lat_diff = float(np.abs(lat[1] - lat[0]))
    expected_first_row = np.round(
        calc_grid_area(
            [
                (lon[0], lat[0] + lat_diff / 2),
                (lon[0], lat[0] - lat_diff / 2),
                (lon[1], lat[0] - lat_diff / 2),
                (lon[1], lat[0] + lat_diff / 2),
            ]
        ),
        2,
    )
    assert np.isclose(area[0, 0], expected_first_row)
