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
import logging
logger = logging.getLogger(__name__)
plt.rcParams["animation.html"] = "jshtml"
from .mask import show

show = show

def ds_reformat_index(ds):
	"""format the dataArray generated from the convert function"""
	if 'lat' in ds.dims and 'lat' in ds.dims:
		return ds.sortby(['lat', 'lon'])
	return ds.reset_coords(['lon', 'lat'], drop = True).rename({'x': 'lon', 'y': 'lat'}).sortby(['lat', 'lon'])


def ds_ts_aggregate(ds, agg_method):
	if agg_method == 'mean':
		ds = ds.transpose('time', 'lat', 'lon').mean(axis = 1).mean(axis = 1)
		return ds
	elif agg_method == 'sum':
		ds = ds.transpose('time', 'lat', 'lon').sum(axis = 1).sum(axis = 1)
		return ds

def time_series(ds, lat_slice = None, lon_slice = None, agg_slice = True,
				agg_slice_method = 'mean', coord_dict = None,
				time_factor = 1, agg_time_method = 'mean',
				figsize = (10, 5), ds_name = None, loc_name = None, title = None, title_size = 12,
				grid = True, legend = True, return_fig = False, **kwargs):
	"""
	Take in the xarray.DataArray, slice of latitude or longitude,
	plot the time series for PM25 at the specified coordinates;
	When users give lat or lon slices, the values can be mean/sum aggregated;
	Users can also provide a dictionary of name-coordinate pairs instead of lat and lon input.
	By default, the method shows the time series aggregated (ds's time unit * time factor)  mean.

	ds (xarray.DataArray): the dataArray object
	lat_slice (tuple): the slice of latitude values
	lon_slice (tuple): the slice of longitude values
	agg_slice (bool): whether the program plot the aggregate values from the slices
	agg_slice_method (str): mean aggregation or sum aggregation for aggregating the spatial slices
	coord_dict (dict): name, coordinate pair of different locations;
			example: {'Beijing': (40, 116.25), 'Shanghai': (31, 121.25)}
	time_factor (float): the factor to aggregate the value of dataArray on
			example: for daily mean on hourly data, time_factor = 24
	agg_time_method (str): mean aggregation or sum aggregation for aggregate time points
	figsize (tuple): the size of the plot
	ds_name (str): name of the DataArray to be shown on title
	loc_name (str): location of the place to be shown on title
	title (str): the title of the result plot
	title_size (float): the size of the title of the result plolt
	grid (bool): add grid lines to the plot, True by default
	legend (bool): add legend to the plot if multiple locations are provided, True by default
	return_fig (bool): return the figure, True by default
	**kwargs: other arguments for xarray.DataArray.plot()
	"""

	ds = ds_reformat_index(ds)
	if not ds_name:
		ds_name = ds.name

	fig = plt.figure(figsize = figsize)
	ax = fig.add_subplot()

	create_title = f"{ds.name} time series - "

	if agg_time_method == 'mean':
		ds = ds.coarsen(time = time_factor, boundary="trim").mean() 
	elif agg_time_method == 'sum':
		ds = ds.coarsen(time = time_factor, boundary="trim").sum() 

	if agg_slice is False and (lat_slice is None and lon_slice is None):
		assert agg_slice == False, "agg_slice cannot be set to False without lat_slice or lon_slice."

	if lat_slice: 
		assert lat_slice[1] > lat_slice[0], "Please give correct latitude slice"
		ds = ds.where(ds.lat >= lat_slice[0], drop = True).where(ds.lat <= lat_slice[1], drop = True)
		create_title += f"lat slice {lat_slice} "

	if lon_slice: 
		assert lon_slice[1] > lon_slice[0], "Please give correct longitude slice"
		ds = ds.where(ds.lon >= lon_slice[0], drop = True).where(ds.lon <= lon_slice[1], drop = True)
		create_title += f"lon slice {lon_slice} "

	if coord_dict:
		if not loc_name:
			loc_name = ', '.join(coord_dict.keys())
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
				logger.info("Find grid cell containing coordinate for %s at lat = %la, lon = %lo.", key, la, lo)
			ds.sel(lat = la, lon = lo).plot(ax = ax, label = key, **kwargs)
			create_title += f"{loc_name} "

	if not coord_dict and (agg_slice is True or (lat_slice is None and lon_slice is None)):
		ds = ds_ts_aggregate(ds, agg_slice_method)
		ds.plot(ax = ax, **kwargs)
		create_title += f"spatially {agg_slice_method} aggregated "

	if agg_slice is False and (lat_slice is not None or lon_slice is not None):
		if lat_slice and not lon_slice:
			if agg_slice_method == 'mean':
				ds = ds.mean(axis = 2)
			elif agg_slice_method == 'sum':
				ds = ds.sum(axis = 2)
			create_title += f"with longitude {agg_slice_method} aggregated "
			for la in ds.lat.values:
				ds.sel(lat = la).plot(ax = ax, label = f"lat {la}", **kwargs)

		elif lon_slice and not lat_slice:
			if agg_slice_method == 'mean':
				ds = ds.mean(axis = 1)
			elif agg_slice_method == 'sum':
				ds = ds.sum(axis = 1)
			create_title += f"with latitude {agg_slice_method} aggregated "
			for lo in ds.lon.values:
				ds.sel(lon = lo).plot(ax = ax, label = f"lon {lo}", **kwargs)

		elif lat_slice and lon_slice:
			for la in ds.lat.values:
				for lo in ds.lon.values:
					ds.sel(lat = la, lon = lo).plot(ax = ax, label = f"lat {la}, lon {lo}", **kwargs)

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

def heatmap(ds, t = None, agg_method = 'mean', shape = None,  shape_width = 0.5, shape_color = 'black',
			map_type = 'colormesh', cmap = 'bone_r', figsize = (10, 6), title = None, title_size = 12,
			grid = True, return_fig = False, **kwargs):
	"""
	Take the xrray dataArray, and a time index or string, plot contour/colormesh map for its values;

	ds: dataArray object
	t: time, either a numeric time index, or the string with the same time format from the dataArray
	m_type: must be either 'contour' or 'colormesh'
	
	ds (xarray.DataArray): the dataArray object
	t (int or str): timepoint to show the map. It can be either a numeric time index,
		or a time string from the xarray.DataArray time dimension
	agg_method (str): if t not provided, perform mean aggregation or sum aggregation
	shape (geopandas.geoseries): shapes to be plotted over the raster.
	shape_width (float): the line width for plotting shapes. 0.5 by default.
	shape_color (str): color of the shape line. Black by default.
	map_type (str): map type, either 'contour' or 'colormesh'
	cmap (str): the color of the heat map, select one from matplotlib.pyplot.colormaps()
	figsize (tuple): the size of the plot
	title (str): the title of the result plot
	title_size (float): the size of the title of the result plolt
	coastlines (bool): add coast lines to the plot, True by default
	grid (bool): add grid lines to the plot, True by default
	return_fig (bool): return the figure, True by default
	**kwargs: other argument for xarray.DataArray.plot.pcolormesh() or
		xarray.DataArray.plot.contourf()
	"""
	
	assert map_type == 'contour' or map_type == 'colormesh', "map_type should either be 'contour' or 'colormesh'"
	assert cmap in plt.colormaps(), ("Please see available colormaps through: matplotlib.pyplot.colormaps() or" +
							   " https://matplotlib.org/stable/gallery/color/colormap_reference.html")

	ds = ds_reformat_index(ds)

	fig = plt.figure(figsize = figsize)
	ax = fig.add_subplot()

	if t is not None:
		if isinstance(t, int):
			time_idx = ds.time.data[t]
			ds = ds.isel(time = t)
		else:
			ds = ds.sel(time = t)
	else:
		if agg_method == 'mean':
			ds = ds.mean(axis = 0)
		elif agg_method == 'sum':
			ds = ds.sum(axis = 0)

	if map_type == 'contour':
		ds.plot.contourf('lon', 'lat', ax=ax, cmap = cmap, **kwargs)
	elif map_type == 'colormesh':
		ds.plot.pcolormesh('lon', 'lat', ax=ax, cmap = cmap, **kwargs)

	if shape is not None:
		shape.boundary.plot(ax=ax, linewidth = shape_width, color = shape_color)

	ax.set_xlim(ds.lon.min(), ds.lon.max())
	ax.set_ylim(ds.lat.min(), ds.lat.max())

	if not title:
		if t == None:
			title = f"{ds.name} aggregated {agg_method}"
		elif isinstance(t, int):
			title = f"{ds.name} Amount at time index {t} - {time_idx}"
		elif isinstance(t, str):
			title = f"{ds.name} Amount at {t}"

	ax.set_title(title, size = title_size)

	if grid:
		ax.grid()

	if return_fig:
		return fig

def heatmap_animation(ds, time_factor = 1, agg_method ='mean',
						shape = None,  shape_width = 0.5, shape_color = 'black',
						cmap = 'bone_r', v_max = None, ds_name = None,
						figsize = (10, 5), title = None, title_size = 12 
						grid = True, **kwargs):
	"""
	Created animated version of colormesh() so users can see the value change over time
	at default, each frame is the average or sum of value per time_unit * time_factor.

	ds (xarray.DataArray): the dataArray object
	time_factor (float): the factor to aggregate the value of dataArray on
			example: for daily mean on hourly data, time_factor = 24
	agg_method (str): mean aggregation or sum aggregation
	shape (geopandas.geoseries): shapes to be plotted over the raster.
	shape_width (float): the line width for plotting shapes. 0.5 by default.
	shape_color (str): color of the shape line. Black by default.
	cmap (str): the color of the heat map, select one from matplotlib.pyplot.colormaps()
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
	assert cmap in plt.colormaps(), ("Please see available colormaps through: matplotlib.pyplot.colormaps() or" +
							   " https://matplotlib.org/stable/gallery/color/colormap_reference.html")

	ds = ds_reformat_index(ds)

	if not ds_name:
		ds_name = ds.name

	if agg_method == 'mean':
		ds = ds.coarsen(time = time_factor, boundary="trim").mean()
	elif agg_method == 'sum':
		ds = ds.coarsen(time = time_factor, boundary="trim").sum()

	if v_max is None: 
		v_max = ds.max()

	fig = plt.figure(figsize= figsize)
	ax = fig.add_subplot()

	if shape is not None:
		shape.boundary.plot(ax=ax, linewidth = shape_width, color = shape_color)

	#initial frame
	image = ds.isel(time = 0).plot.imshow(ax=ax, vmin = 0, vmax = v_max, cmap = cmap, **kwargs)
	if grid:
		ax.grid()

	def update(t):
		"""function to update each frame"""
		if title is None:
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
