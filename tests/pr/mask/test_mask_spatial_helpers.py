import numpy as np
import xarray as xr

from geodata.mask.spatial import calc_grid_area, coarsen, ds_reformat_index


def test_ds_reformat_index_renames_and_sorts_xy():
    x = np.array([101.0, 100.5, 100.0])
    y = np.array([30.0, 30.5, 31.0])
    arr = np.arange(9, dtype=np.float32).reshape(3, 3)
    da = xr.DataArray(arr, dims=("y", "x"), coords={"x": x, "y": y}, name="signal")

    out = ds_reformat_index(da)
    assert out.dims == ("lat", "lon")
    assert np.all(np.diff(out["lat"].values) >= 0)
    assert np.all(np.diff(out["lon"].values) >= 0)


def test_coarsen_mean_on_aligned_grid():
    lat_hi = np.array([0.0, 0.25, 0.5, 0.75])
    lon_hi = np.array([10.0, 10.25, 10.5, 10.75])
    hi = xr.DataArray(
        np.arange(16, dtype=np.float32).reshape(4, 4),
        dims=("lat", "lon"),
        coords={"lat": lat_hi, "lon": lon_hi},
        name="mask",
    )

    lat_lo = np.array([0.125, 0.625])
    lon_lo = np.array([10.125, 10.625])
    lo = xr.Dataset(coords={"lat": lat_lo, "lon": lon_lo})

    out = coarsen(hi, lo, func="mean")
    # Freeze current legacy coarsen behavior.
    expected = np.array([[7.5, 9.0], [13.5, 15.0]], dtype=np.float32)
    np.testing.assert_allclose(out.values, expected, atol=1e-6)


def test_calc_grid_area_positive_and_latitude_sensitive():
    # Avoid perfectly symmetric parallels around 0 that can trip AEA constraints.
    cell_low_lat = [(0.0, 1.5), (0.0, 0.5), (1.0, 0.5), (1.0, 1.5)]
    cell_high_lat = [(0.0, 60.5), (0.0, 59.5), (1.0, 59.5), (1.0, 60.5)]

    area_low_lat = calc_grid_area(cell_low_lat)
    area_high_lat = calc_grid_area(cell_high_lat)

    assert area_low_lat > 0
    assert area_high_lat > 0
    assert area_low_lat > area_high_lat
