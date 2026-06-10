"""Regression tests for merge_layer with in-memory (/vsimem) layers."""

from pathlib import Path

import numpy as np
from rasterio.transform import from_bounds

from geodata.mask import Mask, save_raster


def _write_layer(path: Path, west: float, south: float, east: float, north: float, pattern: str):
    nlon, nlat = 8, 6
    transform = from_bounds(west, south, east, north, nlon, nlat)
    arr = np.zeros((nlat, nlon), dtype=np.uint8)
    if pattern == "left":
        arr[:, : nlon // 2] = 1
    elif pattern == "right":
        arr[:, nlon // 2 :] = 1
    else:
        arr[nlat // 4 : 3 * nlat // 4, nlon // 4 : 3 * nlon // 4] = 1
    save_raster(arr, transform, str(path))
    return transform


def _mask_with_filtered_layers(tmp_path: Path, *, overlap: bool = False) -> Mask:
    west, south, east, north = 100.0, 30.0, 101.0, 31.0
    layer_a = tmp_path / "layer_a.tif"
    layer_b = tmp_path / "layer_b.tif"
    if overlap:
        _write_layer(layer_a, west, south, east, north, "center")
        _write_layer(layer_b, west, south, east, north, "center")
    else:
        _write_layer(layer_a, west, south, east, north, "left")
        _write_layer(layer_b, west, south, east, north, "right")

    mask = Mask("inmemory_merge_test", mask_dir=str(tmp_path / "masks"))
    mask.add_layer(str(layer_a), layer_name="a")
    mask.add_layer(str(layer_b), layer_name="b")
    mask.filter_layer("a", min_bound=0.5, binarize=True, dest_layer_name="a")
    mask.filter_layer("b", min_bound=0.5, binarize=True, dest_layer_name="b")
    return mask


def test_filtered_layers_are_vsimem_backed(tmp_path):
    mask = _mask_with_filtered_layers(tmp_path)
    for ds in mask.layers.values():
        assert ds.name.startswith("/vsimem"), ds.name
        ds.read(1)


def test_merge_and_after_filter_layer(tmp_path):
    mask = _mask_with_filtered_layers(tmp_path)
    merged = mask.merge_layer(
        method="and",
        layers=["a", "b"],
        reference_layer="a",
        show_raster=False,
    )
    assert not merged.closed
    data = merged.read(1)
    assert data.shape == (6, 8)
    assert mask.merged_mask is not None
    assert not mask.saved


def test_merge_sum_after_filter_layer(tmp_path):
    mask = _mask_with_filtered_layers(tmp_path)
    merged = mask.merge_layer(
        method="sum",
        layers=["a", "b"],
        weights={"a": 1.0, "b": 2.0},
        reference_layer="a",
        show_raster=False,
        attribute_save=False,
    )
    assert not merged.closed
    data = merged.read(1)
    assert np.any(data > 0)


def test_merge_and_trim_after_filter(tmp_path):
    mask = _mask_with_filtered_layers(tmp_path, overlap=True)
    merged = mask.merge_layer(
        method="and",
        layers=["a", "b"],
        reference_layer="a",
        trim=True,
        show_raster=False,
    )
    data = merged.read(1)
    assert data.shape[0] <= 6
    assert data.shape[1] <= 8
    assert np.any(data != 0)
