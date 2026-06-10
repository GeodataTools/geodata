from pathlib import Path
from typing import Any, cast

import numpy as np
import pytest
import xarray as xr
from rasterio.transform import from_bounds

from geodata.cutout import Cutout
from geodata.mask import Mask, save_raster


def _build_minimal_cutout() -> Cutout:
    cutout = Cutout.__new__(Cutout)
    cutout.name = "legacy-error-cutout"
    cutout.meta = xr.Dataset(
        coords={
            "x": np.array([100.0, 100.25, 100.5]),
            "y": np.array([30.5, 30.25, 30.0]),
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


def _sample_dataset() -> xr.Dataset:
    t = np.array(["2016-01-01T00:00:00"], dtype="datetime64[ns]")
    y = np.array([30.5, 30.25, 30.0])
    x = np.array([100.0, 100.25, 100.5])
    data = np.arange(len(t) * len(y) * len(x), dtype=np.float32).reshape(
        len(t), len(y), len(x)
    )
    return xr.Dataset(
        {"signal": (("time", "y", "x"), data)},
        coords={"time": t, "y": y, "x": x},
    )


def _create_saved_empty_mask(mask_dir: Path, name: str) -> None:
    # Create a mask object that is saved but has no merged/shape masks.
    mask = Mask(name=name, mask_dir=str(mask_dir))
    mask.save_mask()


def _create_unsaved_mask_with_layer(mask_dir: Path, name: str) -> Mask:
    west, south, east, north = 100.0, 30.0, 100.75, 30.75
    arr = np.ones((3, 3), dtype=np.uint8)
    transform = from_bounds(west, south, east, north, arr.shape[1], arr.shape[0])
    layer_path = mask_dir / f"{name}.tif"
    save_raster(arr, transform, str(layer_path))
    mask = Mask(name=name, mask_dir=str(mask_dir))
    mask.add_layer(str(layer_path), layer_name="base")
    return mask


def test_mask_raises_without_added_masks():
    cutout = _build_minimal_cutout()
    ds = _sample_dataset()

    with pytest.raises(ValueError, match="No mask found in cutout"):
        cutout.mask(ds)


def test_mask_raises_when_true_area_requested_without_area():
    cutout = _build_minimal_cutout()
    ds = _sample_dataset()
    cutout.merged_mask = xr.DataArray(
        np.ones((1, 3, 3), dtype=np.float32),
        dims=("band", "lat", "lon"),
        coords={
            "band": [1],
            "lat": ds["y"].values,
            "lon": ds["x"].values,
        },
    )

    with pytest.raises(ValueError, match="No area data found"):
        cutout.mask(ds, true_area=True)


def test_add_mask_raises_for_saved_mask_without_merged_or_shape(tmp_path, monkeypatch):
    cutout = _build_minimal_cutout()
    mask_dir = tmp_path / "masks"
    mask_dir.mkdir(parents=True, exist_ok=True)
    name = "empty_saved_mask"
    _create_saved_empty_mask(mask_dir, name)

    monkeypatch.setattr("geodata.cutout.config.MASK_DIR", str(mask_dir))

    with pytest.raises(ValueError, match=f"No mask found in {name}"):
        cutout.add_mask(name)


def test_mask_load_xarray_raises_when_unsaved(tmp_path):
    mask_dir = tmp_path / "masks"
    mask_dir.mkdir(parents=True, exist_ok=True)
    mask = _create_unsaved_mask_with_layer(mask_dir, name="unsaved_mask")

    with pytest.raises(ValueError, match="has not been saved"):
        mask.load_merged_xr()

    with pytest.raises(ValueError, match="has not been saved"):
        _ = mask.load_shape_xr(names=cast(Any, []))
