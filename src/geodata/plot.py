# Copyright 2022 Jiahe Feng (Davidson Lab)
# Copyright 2022 Xiqiang Liu (Davidson Lab)

# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License as
# published by the Free Software Foundation; either version 3 of the
# License, or (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

import logging
from typing import Literal, Optional, Union

import geopandas as gpd
import matplotlib.animation as anim
import matplotlib.pyplot as plt
import xarray as xr

from .cutout import ds_reformat_index
from .mask import show  # noqa: F401

plt.rcParams["animation.html"] = "jshtml"

logger = logging.getLogger(__name__)
CoordinateType = tuple[Union[float, int], Union[float, int]]
ReductionType = Literal["mean", "sum"]


def ds_ts_aggregate(
    ds: Union[xr.Dataset, xr.DataArray], agg_method: ReductionType
) -> Union[xr.DataArray, xr.Dataset]:
    """Aggregate the xarray.Dataset or xarray.DataArray along the lat and lon dimensions.

    Args:
        ds (Union[xr.Dataset, xr.DataArray]): The xarray Dataset.
        agg_method (Literal["mean", "sum"]): The aggregation method. If "mean", the mean
            of all of the values will be taken. If "sum", the sum of all of the values will be taken.

    Returns:
        Union[xr.Dataset, xr.DataArray]: The aggregated xarray Dataset.

    Raises:
        NotImplementedError: If the aggregation method is not supported.
    """

    if agg_method == "mean":
        ds = ds.transpose("time", "lat", "lon").mean(axis=1).mean(axis=1)
    elif agg_method == "sum":
        ds = ds.transpose("time", "lat", "lon").sum(axis=1).sum(axis=1)
    else:
        raise NotImplementedError(f"agg_method {agg_method} is not supported.")
    return ds


def time_series(
    ds: xr.DataArray,
    lat_slice: Optional[CoordinateType] = None,
    lon_slice: Optional[CoordinateType] = None,
    agg_slice: bool = True,
    agg_slice_method: ReductionType = "mean",
    coord_dict: Optional[dict[str, CoordinateType]] = None,
    time_factor: float = 1.0,
    agg_time_method: ReductionType = "mean",
    figsize: tuple[float, float] = (10.0, 5.0),
    ds_name: Optional[str] = None,
    loc_name: Optional[str] = None,
    title: Optional[str] = None,
    title_size: float = 12.0,
    grid: bool = True,
    legend: bool = True,
    return_fig: bool = False,
    **kwargs,
) -> Optional[plt.Figure]:
    """Take in the xarray.DataArray, slice of latitude or longitude,
    plot the time series. When users give lat or lon slices, the values can be mean/sum aggregated.
    Users can also provide a dictionary of name-coordinate pairs instead of lat and lon input.
    By default, the method shows the time series' aggregated mean.

    Args:
        ds (xr.DataArray): The target xarray.DataArray object.
        lat_slice (tuple): The slice of latitude values.
        lon_slice (tuple): The slice of longitude values
        agg_slice (bool): Whether the program plot the aggregate values from the slices
        agg_slice_method (str): Reduction method for aggregating the spatial slices. This can be either
            "mean" or "sum".
        coord_dict (dict): The (Name, Coordinate) pair of different locations;
            An example: `{'Beijing': (40, 116.25), 'Shanghai': (31, 121.25)}`
        time_factor (float): The factor to aggregate the value of dataArray on
            An example: for daily mean on hourly data, time_factor = 24
        agg_time_method (str): Reduction method for time dimension. This can be either "mean" or "sum".
        figsize (tuple): The size of the plot
        ds_name (str): Name of the DataArray to be shown on title
        loc_name (str): Location of the place to be shown on title
        title (str): The title of the result plot
        title_size (float): The size of the title of the result plolt
        grid (bool): Whether to add grid lines to the plot, True by default
        legend (bool): Add legend to the plot if multiple locations are provided, True by default
        return_fig (bool): Whether to return the figure or not. True by default
        **kwargs: Other keyword arguments for xarray.DataArray.plot()

    Returns:
        Optional[plt.Figure]: The figure object if `return_fig` is True.

    Raises:
        NotImplementedError: If the reduction method is not supported.
    """

    ds = ds_reformat_index(ds)
    if not ds_name:
        ds_name = ds.name

    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot()

    create_title = f"{ds.name} Time Series - "

    ds = ds.coarsen(time=time_factor, boundary="trim")
    if agg_time_method == "mean":
        ds = ds.mean()
    elif agg_time_method == "sum":
        ds = ds.sum()
    else:
        raise NotImplementedError(
            f"agg_time_method {agg_time_method} is not supported."
        )

    if not agg_slice and lat_slice is None and lon_slice is None:
        if agg_slice:
            raise RuntimeError(
                "agg_slice cannot be set to True without lat_slice or lon_slice."
            )

    if lat_slice:
        if lat_slice[1] < lat_slice[0]:
            raise ValueError(
                "Please give correct latitude slice. The second value should be larger than the first."
            )
        ds = ds.where(ds.lat >= lat_slice[0], drop=True).where(
            ds.lat <= lat_slice[1], drop=True
        )
        create_title += f"lat slice {lat_slice} "

    if lon_slice:
        if lon_slice[1] < lon_slice[0]:
            raise ValueError(
                "Please give correct longitude slice. The second value should be larger than the first."
            )
        ds = ds.where(ds.lon >= lon_slice[0], drop=True).where(
            ds.lon <= lon_slice[1], drop=True
        )
        create_title += f"lon slice {lon_slice} "

    if coord_dict:
        if not loc_name:
            loc_name = ", ".join(coord_dict.keys())
        for key, value in coord_dict.items():
            all_lat = ds.lat.data
            all_lon = ds.lon.data
            la, lo = value[0], value[1]
            log_new_coord = False
            if la not in all_lat:
                if la > all_lat.min() and la < all_lat.max():
                    log_new_coord = True
                    la = all_lat[all_lat < la][-1]
                else:
                    raise ValueError(f"Latitude for {key} out of bound.")

            if lo not in all_lon:
                if lo > all_lon.min() and lo < all_lon.max():
                    log_new_coord = True
                    lo = all_lon[all_lon < lo][-1]
                else:
                    raise ValueError(f"Longitude for {key} out of bound.")

            if log_new_coord:
                logger.info(
                    "Find grid cell containing coordinate for %s at lat = %f, lon = %f.",
                    key,
                    la,
                    lo,
                )
            ds.sel(lat=la, lon=lo).plot(ax=ax, label=key, **kwargs)
            create_title += f"{loc_name} "

    if not coord_dict and (
        agg_slice is True or (lat_slice is None and lon_slice is None)
    ):
        ds = ds_ts_aggregate(ds, agg_slice_method)
        ds.plot(ax=ax, **kwargs)
        create_title += f"spatially {agg_slice_method} aggregated "

    if agg_slice is False and (lat_slice is not None or lon_slice is not None):
        if lat_slice and not lon_slice:
            if agg_slice_method == "mean":
                ds = ds.mean(axis=2)
            elif agg_slice_method == "sum":
                ds = ds.sum(axis=2)
            create_title += f"with longitude {agg_slice_method} aggregated "
            for la in ds.lat.values:
                ds.sel(lat=la).plot(ax=ax, label=f"lat {la}", **kwargs)

        elif lon_slice and not lat_slice:
            if agg_slice_method == "mean":
                ds = ds.mean(axis=1)
            elif agg_slice_method == "sum":
                ds = ds.sum(axis=1)
            create_title += f"with latitude {agg_slice_method} aggregated "
            for lo in ds.lon.values:
                ds.sel(lon=lo).plot(ax=ax, label=f"lon {lo}", **kwargs)

        elif lat_slice and lon_slice:
            for la in ds.lat.values:
                for lo in ds.lon.values:
                    ds.sel(lat=la, lon=lo).plot(
                        ax=ax, label=f"lat {la}, lon {lo}", **kwargs
                    )

    if legend and (agg_slice is False or coord_dict):
        ax.legend()
    if grid:
        ax.grid()

    if time_factor > 1:
        create_title += f"- time aggregated by factor of {time_factor}."
    if not title:
        title = create_title
    ax.set_title(title, size=title_size)

    if return_fig:
        return fig


def heatmap(
    ds: xr.DataArray,
    t: Optional[Union[int, str]] = None,
    agg_method: ReductionType = "mean",
    shape: Optional[gpd.GeoSeries] = None,
    shape_width: float = 0.5,
    shape_color: str = "black",
    map_type: Literal["contour", "colormesh"] = "colormesh",
    cmap: str = "bone_r",
    figsize: tuple[float, float] = (10.0, 6.0),
    title: Optional[str] = None,
    title_size: float = 12,
    grid: bool = True,
    return_fig: bool = False,
    **kwargs,
) -> Optional[plt.Figure]:
    """Take an xarray.DataArray and a time index or string, plot contour/colormesh map for its values.

    Args:
        ds (xr.DataArray): The target DataArray object.
        t (Optional[Union[int, str]]): Target timestamp. This could either a numeric time index,
            or a time string from the xarray.DataArray time dimension.
        agg_method (Literal["mean", "sum"]): Aggregation method in the time dimension.
            This is used if t was not not provided. Options can either be mean aggregation or sum aggregation.
        shape (geopandas.GeoSeries): Shapes to be plotted over the raster.
        shape_width (float): Width of lines for plotting shapes. 0.5 by default.
        shape_color (str): Color of the shape line. Black by default.
        map_type (Literal["contour", "colormesh"]): Map type. This can either be "contour" or "colormesh".
        cmap (str): The color of the heat map, select one from matplotlib.pyplot.colormaps.
        figsize (tuple): The size of the plot.
        title (Optional[str]): The title of the result plot. Optional. If not provided, the title will be
            automatically generated.
        title_size (float): The size of the title of the result plot. 12 by default.
        coastlines (bool): Whether to add coast lines to the plot, True by default.
        grid (bool): Whether to add grid lines to the plot, True by default.
        return_fig (bool): Whether to return the plt.Figure object. True by default
        **kwargs: Additional arguments for xarray.DataArray.plot.pcolormesh or
            xarray.DataArray.plot.contourf, depending on selected `map_type`.

    Returns:
        Optional[plt.Figure]: The figure object if `return_fig` is True.

    Raises:
        ValueError: If the map type is not supported.
    """

    if map_type not in {"contour", "colormesh"}:
        raise ValueError(f"map_type {map_type} is not supported.")
    if cmap not in plt.colormaps():
        raise ValueError(
            "Please see available colormaps through: matplotlib.pyplot.colormaps() or "
            "https://matplotlib.org/stable/gallery/color/colormap_reference.html"
        )

    ds = ds_reformat_index(ds)

    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot()

    if t is not None:
        if isinstance(t, int):
            time_idx = ds.time.data[t]
            ds = ds.isel(time=t)
        else:
            ds = ds.sel(time=t)
    else:
        if agg_method == "mean":
            ds = ds.mean(axis=0)
        elif agg_method == "sum":
            ds = ds.sum(axis=0)

    if map_type == "contour":
        ds.plot.contourf("lon", "lat", ax=ax, cmap=cmap, **kwargs)
    elif map_type == "colormesh":
        ds.plot.pcolormesh("lon", "lat", ax=ax, cmap=cmap, **kwargs)

    if shape is not None:
        shape.boundary.plot(ax=ax, linewidth=shape_width, color=shape_color)

    ax.set_xlim(ds.lon.min(), ds.lon.max())
    ax.set_ylim(ds.lat.min(), ds.lat.max())

    if not title:
        if t is None:
            title = f"{ds.name} aggregated {agg_method}"
        elif isinstance(t, int):
            title = f"{ds.name} Amount at time index {t} - {time_idx}"
        elif isinstance(t, str):
            title = f"{ds.name} Amount at {t}"

    ax.set_title(title, size=title_size)

    if grid:
        ax.grid()

    if return_fig:
        return fig


def heatmap_animation(
    ds: xr.DataArray,
    time_factor: float = 1,
    agg_method: ReductionType = "mean",
    shape: Optional[gpd.GeoSeries] = None,
    shape_width: float = 0.5,
    shape_color: str = "black",
    cmap: str = "bone_r",
    v_max: Optional[float] = None,
    ds_name: Optional[str] = None,
    figsize: tuple[float, float] = (10, 5),
    title: Optional[str] = None,
    title_size: float = 12,
    grid: bool = True,
    **kwargs,
):
    """Created animated version of `colormesh` so users can see the value change over time
    at default, each frame is the average or sum of value per time_unit * time_factor.

    Args:
        ds (xarray.DataArray): The target DataArray object.
        time_factor (float): Tthe factor to aggregate the value of DataArray on
            Example: for daily mean on hourly data, time_factor = 24. Defaults to 1.
        agg_method (str): Aggregation method. Can be either `mean` or `sum`.
        shape (geopandas.GeoSeries): Shapes to be plotted over the raster.
        shape_width (float): The line width for plotting shapes. 0.5 by default.
        shape_color (str): Color of the shape line. Black by default.
        cmap (str): The color of the heat map, select one from matplotlib.pyplot.colormaps.
        v_max (float): The maximum value in the heatmap.
        ds_name (str): Name of the DataArray to be shown on title.
        figsize (tuple): The size of the plot. (10, 5) by default.
        title (str): The title of the result plot. Optional. If not provided, the title will be
            automatically generated.
        title_size (float): The size of the title of the result plot.
        coastlines (bool): Whether to add coast lines to the plot, True by default.
        grid (bool): Whether to add grid lines to the plot, True by default.
        **kwargs (dict): Additional arguments for xarray.DataArray.plot.imshow.

    Returns:
        matplotlib.animation.FuncAnimation: The animation object.

    Raises:
        ValueError: If the DataArray does not contain the time dimension.
        ValueError: If the colormap is not supported.
        NotImplementedError: If the aggregation method is not supported.
    """

    if "time" not in ds.dims:
        raise ValueError("The DataArray must contain the time dimension")
    if cmap not in plt.colormaps():
        raise ValueError(
            "Please see available colormaps through: matplotlib.pyplot.colormaps() or "
            "https://matplotlib.org/stable/gallery/color/colormap_reference.html"
        )

    ds = ds_reformat_index(ds)

    if not ds_name:
        ds_name = ds.name

    if agg_method == "mean":
        ds = ds.coarsen(time=time_factor, boundary="trim").mean()
    elif agg_method == "sum":
        ds = ds.coarsen(time=time_factor, boundary="trim").sum()
    else:
        raise NotImplementedError(f"agg_method {agg_method} is not supported.")

    if v_max is None:
        v_max = ds.max()

    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot()

    if shape is not None:
        shape.boundary.plot(ax=ax, linewidth=shape_width, color=shape_color)

    # initial frame
    image = ds.isel(time=0).plot.imshow(ax=ax, vmin=0, vmax=v_max, cmap=cmap, **kwargs)
    if grid:
        ax.grid()

    def update(t):
        """function to update each frame"""
        if title is None:
            ax.set_title(f"{ds_name} at time = {t}", size=title_size)
        else:
            ax.set_title(title, size=title_size)
        image.set_array(ds.sel(time=t))
        return image

    animation = anim.FuncAnimation(fig, update, frames=ds.time.values, blit=False)

    return animation


def save_animation(file_name: str):
    """If the Ipython notebook is opened in a browser, and an animation output was already generated.
    This functioon saves the animation to a file.

    Args:
        file_name (str): The output file name.
    """

    javascript = (
        """
    <script type="text/Javascript">
        function set_value(){
            elements = document.getElementsByClassName('output_subarea output_html rendered_html output_result')
            var var_values = ''
            for (i = 0; i < elements.length; i++){
                if (elements[i].getElementsByClassName('animation').length != 0){
                var_values += elements[i].innerHTML;
            }}

            (function(console){
            /* credit of the console.save function: stackoverflow.com/questions/11849562/*/
            console.save = function(data, filename){
                if(!data) {
                    console.error('Console.save: No data')
                    return;
                }
                if(!filename) filename = 'console.json'
                if(typeof data === "object"){
                    data = JSON.stringify(data, undefined, 4)
                }
                var blob = new Blob([data], {type: 'text/json'}),
                    e    = document.createEvent('MouseEvents'),
                    a    = document.createElement('a')
                a.download = filename
                a.href = window.URL.createObjectURL(blob)
                a.dataset.downloadurl =  ['text/json', a.download, a.href].join(':')
                e.initMouseEvent('click', true, false, window, 0, 0, 0, 0, 0, false, false, false, false, 0, null)
                a.dispatchEvent(e)
            }
            })(console)
            console.save(var_values, '"""  # noqa: E501
        + file_name
        + """')
        }
        set_value()
    </script>
    """
    )
    return javascript
