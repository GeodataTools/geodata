"""Shared spatial helper utilities for masking workflows."""

from functools import partial
from typing import Any, Literal, cast

import numpy as np
import pyproj
import shapely
import xarray as xr
from shapely import ops


def ds_reformat_index(ds: xr.Dataset | xr.DataArray) -> xr.Dataset | xr.DataArray:
    """Normalize data coordinates to sorted ``lat``/``lon``."""
    if "lat" in ds.dims and "lon" in ds.dims:
        return ds.sortby(["lat", "lon"])
    if "lat" in ds.coords and "lon" in ds.coords:
        return (
            ds.reset_coords(["lon", "lat"], drop=True)
            .rename({"x": "lon", "y": "lat"})
            .sortby(["lat", "lon"])
        )
    return ds.rename({"x": "lon", "y": "lat"}).sortby(["lat", "lon"])


def _find_intercept(list1, list2, start, threshold=0):
    """Find best start offset for coarsening alignment."""
    min_res = 0
    init = 0
    i = 0
    for i in range(len(list1) - start):
        resid = ((list1[start + i] - list2[0]) % (list2[1] - list2[0])).values.tolist()
        if i == 0:
            init = resid
        if resid <= threshold:
            return i
        if resid > min_res:
            min_res = resid
        else:
            min_res = resid
            break
    if min_res == init:
        return 0
    return i


def coarsen(
    ori: xr.Dataset | xr.DataArray,
    tar: xr.Dataset | xr.DataArray,
    func: Literal["sum", "mean"] = "mean",
):
    """Reindex/coarsen ``ori`` according to target coordinates in ``tar``."""
    lat_multiple = round(
        ((tar.lat[1] - tar.lat[0]) / (ori.lat[1] - ori.lat[0])).values.tolist()
    )
    lon_multiple = round(
        ((tar.lon[1] - tar.lon[0]) / (ori.lon[1] - ori.lon[0])).values.tolist()
    )
    lat_start = _find_intercept(ori.lat, tar.lat, (lat_multiple - 1) // 2)
    lon_start = _find_intercept(ori.lon, tar.lon, (lon_multiple - 1) // 2)

    if func == "mean":
        coarsened = (
            ori.isel(lat=slice(lat_start, None), lon=slice(lon_start, None))
            .coarsen(
                dim={"lat": lat_multiple, "lon": lon_multiple},
                side={"lat": "left", "lon": "left"},
                boundary="pad",
            )
        )
        reduced = cast(Any, coarsened).mean()
    elif func == "sum":
        coarsened = (
            ori.isel(lat=slice(lat_start, None), lon=slice(lon_start, None))
            .coarsen(
                dim={"lat": lat_multiple, "lon": lon_multiple},
                side={"lat": "left", "lon": "left"},
                boundary="pad",
            )
        )
        reduced = cast(Any, coarsened).sum()
    else:
        raise ValueError("func can only be 'mean' or 'sum'")

    return reduced.reindex_like(tar, method="nearest")


def calc_grid_area(lis_lats_lons):
    """Calculate area in km^2 for a grid cell defined by corner coordinates."""
    lons, lats = zip(*lis_lats_lons)
    ll = list(set(lats))[::-1]
    var = []
    for i in range(len(ll)):
        var.append("lat_" + str(i + 1))
    st = ""
    for v, l in zip(var, ll):  # noqa: E741
        st = st + str(v) + "=" + str(l) + " " + "+"
    st = (
        st
        + "lat_0="
        + str(np.mean(ll))
        + " "
        + "+"
        + "lon_0"
        + "="
        + str(np.mean(lons))
    )
    tx = "+proj=aea +" + st
    pa = pyproj.Proj(tx)

    x, y = pa(lons, lats)
    cop = {"type": "Polygon", "coordinates": [zip(x, y)]}
    return shapely.geometry.shape(cop).area / 1000000


def calc_shp_area(shp, shp_projection="+proj=latlon"):
    """Calculate area in km^2 for a shape object."""
    temp_shape = ops.transform(
        partial(
            pyproj.transform,
            pyproj.Proj(shp_projection),
            pyproj.Proj(proj="aea", lat_1=shp.bounds[1], lat_2=shp.bounds[3]),
        ),
        shp,
    )
    return temp_shape.area / 1000000


__all__ = ["ds_reformat_index", "coarsen", "calc_grid_area", "calc_shp_area"]
