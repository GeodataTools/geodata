## Copyright 2022 Jiahe Feng (Davidson Lab).

## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License as
## published by the Free Software Foundation; either version 3 of the
## License, or (at your option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.

## You should have received a copy of the GNU General Public License
## along with this program. If not, see <http://www.gnu.org/licenses/>.

import matplotlib.animation as anim
import matplotlib.pyplot as plt
plt.rcParams["animation.html"] = "jshtml"
from .mask import show

#packages that requires users to install externally
import cartopy.crs as ccrs

show = show

def ds_reformat_index(ds):
    """format the dataArray generated from the convert function"""
    return ds.reset_coords(['lon', 'lat'], drop = True).rename({'x': 'lon', 'y': 'lat'}).sortby(['lat', 'lon'])


def time_series(ds, lat = None, lon = None, coord_dict = None, time_factor = 1, agg_method = 'mean', 
                figsize = (10, 5), ds_name = None, loc_name = None, title = None, title_size = 12, 
                grid = True, legend = True, return_fig = False, **kwargs):
    """
    Take in the xarray.DataArray, latitude and longitudes, plot the time series for PM25 at the specified coordinates;
    Lat and lon can both be a number, or list, but they can't both be list.
    Users can also provide a dictionary of name-coordinate pairs instead of lat and lon input.
    By default, the method shows the time series aggregated (ds's time unit * time factor)  mean.
    
    ds (xarray.DataArray): the dataArray object
    lat (int or list): a number of a list of latitude
    lon (int or list): a number of a list of longitudes
    coord_dict (dict): name, coordinate pair of different locations; 
            example: {'Beijing': (40, 116.25), 'Shanghai': (31, 121.25)}
    time_factor (float): the factor to aggregate the value of dataArray on
            example: for daily mean on hourly data, time_factor = 24
    agg_method (str): mean aggregation or sum aggregation
    figsize (tuple): the size of the plot
    ds_name (str): name of the DataArray to be shown on title
    loc_name (str): location of the place to be shown on title
    title (str): the title of the result plot
    title_size (float): the size of the title of the result plolt
    grid (bool): add grid lines to the plot, True by default
    legend (bool): add legend to the plot if multiple locations are provided, True by default
    return_fig (bool): return the figure, True by default
    **kwargs (dict): other argument for xarray.DataArray.plot()
    """

    ds = ds_reformat_index(ds)
    if not ds_name: ds_name = ds.name

    fig = plt.figure(figsize = figsize)
    ax = fig.add_subplot()

    if agg_method == 'mean':
        ds = ds.coarsen(time = time_factor, boundary="trim").mean() 
    elif agg_method == 'sum':
        ds = ds.coarsen(time = time_factor, boundary="trim").sum() 

    if lon is None and lat is None and coord_dict is None:
        if agg_method == 'mean':
            ds = ds.transpose('time', 'lat', 'lon').mean(axis = 1).mean(axis = 1)
        elif agg_method == 'sum':
            ds = ds.transpose('time', 'lat', 'lon').sum(axis = 1).sum(axis = 1)

        ds.plot(ax = ax, **kwargs)

    else:
        if coord_dict == None:
            assert (lat != None and lon != None)
            assert not (isinstance(lat, list) and isinstance(lon, list))
        else:
            assert (lat == None and lon == None)

        if coord_dict:
            if not loc_name: loc_name = ', '.join(coord_dict.keys())
            for key, value in coord_dict.items():
                ds.sel(lat = value[0], lon = value[1]).plot(ax = ax, label = key, **kwargs)
            if not title: 
                title = '{} time series at {}'.format(ds_name, loc_name, lat, lon)
        else:
            if not loc_name: 
                loc_name = 'location'
            ds.sel(lat = lat, lon = lon).plot.line(ax = ax, x="time", **kwargs)
            if not title: 
                title = '{} time series at {} - lat {}, lon {}'.format(ds_name, loc_name, lat, lon)
    
    if legend and coord_dict: ax.legend()
    if grid: ax.grid()
        
    ax.set_title(title, size=title_size)
    
    if return_fig: return fig

def heatmap(ds, t, map_type = 'colormesh', color = None, figsize = (10, 6), title = None, title_size = 12, 
                coastlines = True, grid = True, return_fig = False, **kwargs):
    """
    Take the xrray dataArray, and a time index or string, plot contour/colormesh map for its values; 

    ds: dataArray object
    t: time, either a numeric time index, or the string with the same time format from the dataArray
    m_type: must be either 'contour' or 'colormesh'
    
    ds (xarray.DataArray): the dataArray object
    t (int or str): timepoint to show the map. It can be either a numeric time index, 
        or a time string from the xarray.DataArray time dimension
    map_type (str): map type, either 'contour' or 'colormesh'
    color (str): the color of the heat map, select one from matplotlib.pyplot.colormaps()
    figsize (tuple): the size of the plot
    title (str): the title of the result plot
    title_size (float): the size of the title of the result plolt
    coastlines (bool): add coast lines to the plot, True by default
    grid (bool): add grid lines to the plot, True by default
    return_fig (bool): return the figure, True by default
    **kwargs (dict): other argument for xarray.DataArray.plot.pcolormesh() or 
        xarray.DataArray.plot.contourf()
    """
    
    assert map_type == 'contour' or map_type == 'colormesh', "map_type should either be 'contour' or 'colormesh'"
    assert color is not None, ("Please see available colormaps through: matplotlib.pyplot.colormaps() or" + 
                               " https://matplotlib.org/stable/gallery/color/colormap_reference.html")

    fig = plt.figure(figsize = figsize)
    if map_type == 'contour':
        ax = fig.add_subplot(projection=ccrs.Miller(11625))
        ax.set_extent([ds.lon.min(), ds.lon.max(), ds.lat.min(), ds.lat.max()], ccrs.PlateCarree())
        if isinstance(t, int):
            ds.isel(time = t).plot.contourf(ax=ax, transform=ccrs.PlateCarree(), cmap = color, **kwargs)
        else:
            ds.sel(time = t).plot.contourf(ax=ax, transform=ccrs.PlateCarree(), cmap = color, **kwargs)
    elif map_type == 'colormesh':
        ax = fig.add_subplot(projection=ccrs.PlateCarree())
        ax.set_extent([ds.lon.min(), ds.lon.max(), ds.lat.min(), ds.lat.max()], ccrs.PlateCarree())
        if isinstance(t, int):
            ds.isel(time = t).plot.pcolormesh('lon', 'lat', ax=ax, cmap = color, **kwargs)
        else:
            ds.sel(time = t).plot.pcolormesh('lon', 'lat', ax=ax, cmap = color, **kwargs)

    if coastlines:
        ax.coastlines() 
    if not title:
        title = "{} Amount at time: {}".format(ds.name, t)
    ax.set_title(title, size = title_size)
    if grid: 
        ax.gridlines(draw_labels=True)
    
    if return_fig:
        return fig

def heatmap_animation(ds, time_factor = 1, agg_method ='mean', color = None, v_max = None, ds_name = None, 
                        figsize = (10, 5), title = None, title_size = 12, 
                        coastlines = True, grid = True, **kwargs):
    """
    Created animated version of colormesh() so users can see the value change over time
    at default, each frame is the average or sum of value per time_unit * time_factor.
    
    ds (xarray.DataArray): the dataArray object
    time_factor (float): the factor to aggregate the value of dataArray on
            example: for daily mean on hourly data, time_factor = 24
    agg_method (str): mean aggregation or sum aggregation
    color (str): the color of the heat map, select one from matplotlib.pyplot.colormaps()
    v_max (float): the maximum value in the heatmap
    ds_name (str): name of the DataArray to be shown on title
    figsize (tuple): the size of the plot
    title (str): the title of the result plot
    title_size (float): the size of the title of the result plolt
    coastlines (bool): add coast lines to the plot, True by default
    grid (bool): add grid lines to the plot, True by default
    **kwargs (dict): other argument for xarray.DataArray.plot.imshow()
      
    """
    assert 'time' in ds.dims, "The dataArray must contain the time dimension"
    assert color is not None, ("Please see available colormaps through: matplotlib.pyplot.colormaps() or" + 
                               " https://matplotlib.org/stable/gallery/color/colormap_reference.html")
    
    ds = ds_reformat_index(ds)

    if not ds_name: 
        ds_name = ds.name
        
    if agg_method == 'mean':
        ds = ds.coarsen(time = time_factor, boundary="trim").mean() 
    elif agg_method == 'sum':
        ds = ds.coarsen(time = time_factor, boundary="trim").sum() 
    
    if v_max == None: v_max = ds.max()
    
    fig = plt.figure(figsize= figsize)
    ax = plt.axes(projection=ccrs.PlateCarree())
    
    #initial frame
    image = ds.isel(time = 0).plot.imshow(ax=ax, vmin = 0, vmax = v_max, 
                                          transform=ccrs.PlateCarree(), cmap = color,
                                          **kwargs)
    if coastlines:
        ax.coastlines()
    if grid:
        ax.gridlines(draw_labels=True)
    
    def update(t):
        """function to update each frame"""
        if title == None:
            ax.set_title(f"{ds_name} at time = {t}", size = title_size)
        else:
            ax.set_title(title, size = title_size)
        image.set_array(ds.sel(time=t))
        return image

    animation = anim.FuncAnimation(fig, update, frames=ds.time.values, blit=False)

    return animation

def save_animation(file_name):
    '''
    If the Ipython notebook is opened in a browser, and an animation output was already generated;
    Save the animation to a file.
    
    file_name (str): the output file name
    '''
    javascript = """
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
            console.save(var_values, '""" + file_name + """')
        }
        set_value()
    </script>
    """
    return javascript

