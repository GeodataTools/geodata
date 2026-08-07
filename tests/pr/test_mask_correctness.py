"""Offline correctness tests for ``src/geodata/mask.py``.

These tests exercise the mask-module fixes with small synthetic GeoTIFFs
(in-memory or in temporary directories). No network access or downloaded
datasets are required -- only rasterio, rioxarray, numpy, and their
dependencies.

The file is written pytest-style but is also directly executable without
pytest:

    python tests/pr/test_mask_correctness.py
"""

import os
import sys
import tempfile
from pathlib import Path

# Make the in-repo sources importable and keep geodata's storage root out of
# the user's home directory. Both must happen BEFORE geodata is imported.
_REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(_REPO_ROOT / "src"))
os.environ.setdefault("GEODATA_ROOT", tempfile.mkdtemp(prefix="geodata_test_root_"))
os.environ.setdefault("MPLBACKEND", "Agg")

import numpy as np  # noqa: E402
import rasterio as ras  # noqa: E402
from rasterio.transform import from_origin  # noqa: E402

from geodata.mask import (  # noqa: E402
    Mask,
    create_temp_tif,
    filter_raster,
    trim_raster,
)

_TRANSFORM = from_origin(100.0, 10.0, 0.1, 0.1)


def _write_tif(path, arr, nodata=None, crs="EPSG:4326", transform=_TRANSFORM):
    """Write a small single-band GeoTIFF and return its path."""
    arr = np.asarray(arr)
    with ras.open(
        str(path),
        "w",
        driver="GTiff",
        height=arr.shape[0],
        width=arr.shape[1],
        count=1,
        dtype=arr.dtype,
        crs=crs,
        transform=transform,
        nodata=nodata,
    ) as dst:
        dst.write(arr, 1)
    return str(path)


def test_filter_raster_min_and_max_bound_combined():
    """min_bound and max_bound together must AND against the ORIGINAL data.

    Elevation vector [100, 250, 5000, 3000, 8000, 50] with min_bound=200 and
    max_bound=4000 must yield [0, 1, 0, 1, 0, 0]; the historical boolean
    chaining bug returned [1, 1, 1, 1, 1, 1].
    """
    elevation = np.array([[100.0, 250.0, 5000.0, 3000.0, 8000.0, 50.0]])
    raster = create_temp_tif(elevation, _TRANSFORM)

    binarized = filter_raster(raster, min_bound=200, max_bound=4000, binarize=True)
    assert np.array_equal(
        binarized.read(1).flatten(), np.array([0, 1, 0, 1, 0, 0])
    ), f"binarized filter returned {binarized.read(1).flatten().tolist()}"

    # binarize=False must keep the original values of the passing pixels.
    unbinarized = filter_raster(raster, min_bound=200, max_bound=4000, binarize=False)
    assert np.array_equal(
        unbinarized.read(1).flatten(),
        np.array([0.0, 250.0, 0.0, 3000.0, 0.0, 0.0]),
    ), f"unbinarized filter returned {unbinarized.read(1).flatten().tolist()}"


def test_trim_raster_keeps_last_data_row_and_column():
    """The docstring example must trim to shape (5, 3), keeping the last
    nonzero row (row 4) and column (column 4).

    The historical off-by-one fed inclusive stop indices into end-exclusive
    Window slices, cropping away the last data row/column -> shape (4, 2).
    """
    arr = np.array(
        [
            [0, 0, 9, 0, 0, 0, 0],
            [0, 0, 1, 2, 0, 0, 0],
            [0, 0, 2, 3, 4, 0, 0],
            [0, 0, 0, 0, 0, 0, 0],
            [0, 0, 1, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0],
        ],
        dtype=np.uint8,
    )
    raster = create_temp_tif(arr, _TRANSFORM)
    trimmed = trim_raster(raster).read(1)

    assert trimmed.shape == (5, 3), f"trimmed shape is {trimmed.shape}, expected (5, 3)"
    assert np.array_equal(trimmed, arr[0:5, 2:5]), f"trimmed content is\n{trimmed}"
    # The '4' (row 2, col 4) and the row-4 '1' must survive the trim.
    assert 4 in trimmed
    assert trimmed[4].any()


def test_merge_layer_and_all_zero_layer_stays_zero():
    """AND with a fully-excluded (all-zero) layer must return all zeros.

    The historical len(np.unique)==1 first-write heuristic fired on the
    uniform zero window and copied the next layer wholesale, resurrecting
    excluded land (result all ones).
    """
    tmp = tempfile.mkdtemp(prefix="geodata_mask_and_")
    zeros = np.zeros((4, 4), dtype=np.uint8)
    ones = np.ones((4, 4), dtype=np.uint8)

    m = Mask(
        "test_and",
        layer_path={
            "excluded": _write_tif(Path(tmp) / "excluded.tif", zeros),
            "available": _write_tif(Path(tmp) / "available.tif", ones),
        },
        mask_dir=os.path.join(tmp, "masks"),
    )
    merged = m.merge_layer(method="and", show_raster=False).read(1)
    assert (merged == 0).all(), f"AND resurrected excluded pixels:\n{merged}"

    # Control: AND of two mixed binary layers is the elementwise AND.
    tmp2 = tempfile.mkdtemp(prefix="geodata_mask_and2_")
    a = np.array([[1, 1], [0, 0]], dtype=np.uint8)
    b = np.array([[1, 0], [1, 0]], dtype=np.uint8)
    m2 = Mask(
        "test_and_mixed",
        layer_path={
            "a": _write_tif(Path(tmp2) / "a.tif", a),
            "b": _write_tif(Path(tmp2) / "b.tif", b),
        },
        mask_dir=os.path.join(tmp2, "masks"),
    )
    merged2 = m2.merge_layer(method="and", show_raster=False).read(1)
    assert np.array_equal(
        merged2, np.array([[1, 0], [0, 0]])
    ), f"AND of mixed layers returned\n{merged2}"


def test_merge_layer_sum_of_two_ones_is_two():
    """SUM of two all-one layers (default weights of 1) must be 2 everywhere.

    The historical heuristic treated the uniform already-merged window as a
    first write and discarded the accumulated values, returning 1.
    """
    tmp = tempfile.mkdtemp(prefix="geodata_mask_sum_")
    ones = np.ones((4, 4), dtype=np.uint8)

    m = Mask(
        "test_sum",
        layer_path={
            "first": _write_tif(Path(tmp) / "first.tif", ones),
            "second": _write_tif(Path(tmp) / "second.tif", ones),
        },
        mask_dir=os.path.join(tmp, "masks"),
    )
    merged = m.merge_layer(method="sum", show_raster=False).read(1)
    assert (merged == 2).all(), f"SUM of two all-one layers returned\n{merged}"


def test_add_layer_leaves_source_file_nodata_untouched():
    """add_layer must never rewrite the source GeoTIFF on disk.

    The historical implementation opened the user's file in "r+" mode and
    persisted ``nodata = 0`` into it, re-labeling real fills (e.g. -9999) as
    valid data for every other tool. The layer geodata holds in memory must
    instead have the original-nodata pixels masked out to 0.
    """
    tmp = tempfile.mkdtemp(prefix="geodata_mask_nodata_")
    slope = np.array([[1.0, 2.0], [-9999.0, 3.0]], dtype=np.float64)
    path = _write_tif(Path(tmp) / "slope.tif", slope, nodata=-9999.0)

    m = Mask("test_nodata", layer_path=path, layer_name="slope", mask_dir=os.path.join(tmp, "masks"))

    with ras.open(path) as src:
        assert src.nodata == -9999.0, f"source nodata was rewritten to {src.nodata}"
        assert np.array_equal(src.read(1), slope), "source data was rewritten"

    # In-memory normalization: the original nodata pixel is masked to 0 so it
    # can never pass a filter as valid data.
    layer = m.layers["slope"].read(1)
    assert layer[1, 0] == 0, f"original nodata pixel leaked into the layer as {layer[1, 0]}"
    assert layer[0, 0] == 1.0 and layer[1, 1] == 3.0


def test_remerge_after_save_serves_new_merge():
    """merge_layer must invalidate the saved flag: after a re-merge,
    save_mask() + load_merged_xr() must return the NEW merge.

    Historically merge_layer never reset ``saved``, so save_mask() was a
    silent no-op and every downstream consumer got the previous merge from
    disk.
    """
    tmp = tempfile.mkdtemp(prefix="geodata_mask_stale_")
    mask_dir = os.path.join(tmp, "masks")
    # float64 so the fractional re-merge weight cannot be truncated by an
    # integer merge destination (merge adopts the reference layer's dtype).
    ones = np.ones((4, 4), dtype=np.float64)

    m = Mask(
        "test_stale",
        layer_path={
            "a": _write_tif(Path(tmp) / "a.tif", ones),
            "b": _write_tif(Path(tmp) / "b.tif", ones),
        },
        mask_dir=mask_dir,
    )
    m.merge_layer(method="sum", show_raster=False)  # 1 + 1 = 2 everywhere
    m.save_mask()

    # Reload from disk (sets saved=True), then re-merge with a new weight.
    m2 = Mask.from_name("test_stale", mask_dir=mask_dir)
    m2.merge_layer(method="sum", weights={"a": 0.5}, show_raster=False)  # 0.5 + 1 = 1.5
    m2.save_mask()

    loaded = m2.load_merged_xr()
    assert np.allclose(
        loaded.values, 1.5
    ), f"stale merge served: expected 1.5 everywhere, got values {np.unique(loaded.values)}"


if __name__ == "__main__":
    import traceback

    all_tests = [
        (name, fn)
        for name, fn in sorted(globals().items())
        if name.startswith("test_") and callable(fn)
    ]
    failed = []
    for test_name, test_fn in all_tests:
        try:
            test_fn()
            print(f"PASS {test_name}")
        except Exception:
            failed.append(test_name)
            print(f"FAIL {test_name}")
            traceback.print_exc()
    print(f"\n{len(all_tests) - len(failed)}/{len(all_tests)} tests passed")
    sys.exit(1 if failed else 0)
