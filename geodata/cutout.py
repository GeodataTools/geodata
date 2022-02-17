## Copyright 2020 Michael Davidson (UCSD), William Honaker, Jiahe Feng (UCSD), Yuanbo Shi. 
## Copyright 2016-2017 Gorm Andresen (Aarhus University), Jonas Hoersch (FIAS), Tom Brown (FIAS)

## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License as
## published by the Free Software Foundation; either version 3 of the
## License, or (at your option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.

## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.


"""
GEODATA

Geospatial Data Collection and "Pre-Analysis" Tools
"""

from __future__ import absolute_import

import xarray as xr
import numpy as np
import os, sys
from six import iteritems
import pyproj
import shapely
from functools import partial

import logging
logger = logging.getLogger(__name__)

from . import config, datasets

from .convert import (convert_cutout, heat_demand, temperature,
					  wind, pv, solar_thermal, soil_temperature)
from .preparation import (cutout_do_task, cutout_prepare,
						  cutout_produce_specific_dataseries,
						  cutout_get_meta, cutout_get_meta_view)
from .mask import load_mask

class Cutout(object):
	def __init__(self, name=None, cutout_dir=config.cutout_dir, **cutoutparams):
		self.name = name
		self.cutout_dir = os.path.join(cutout_dir, name)
		self.prepared = False
		self.empty = False
		self.meta_append = 0
		self.config = cutoutparams.pop('weather_data_config')
		self.meta = meta = None
		self.merged_mask = None
		self.shape_mask = None
		self.area = None


		if 'bounds' in cutoutparams:
			# if passed bounds array instead of xs, ys slices
			x1, y1, x2, y2 = cutoutparams.pop('bounds')
			cutoutparams.update(xs=slice(x1, x2),
								ys=slice(y1, y2))

		if not 'years' in cutoutparams:
			raise ValueError("`years` need to be specified")
		if not isinstance(cutoutparams['years'], slice):
			years = cutoutparams.pop('years')
			if isinstance(years, list):
				cutoutparams.update(years=slice(years[0], years[-1]))
			else:
				raise TypeError('Unrecognized years parameter given. Only slice or array accepted.')

		if not 'months' in cutoutparams:
			logger.info("No months specified, defaulting to 1-12")
			cutoutparams.update(months=slice(1, 12))

		if os.path.isdir(self.cutout_dir):
			# If cutout dir exists, check completness of files
			if os.path.isfile(self.datasetfn()):  # open existing meta file
				self.meta = meta = xr.open_dataset(self.datasetfn()).stack(**{'year-month': ('year', 'month')})

			if not meta is None and 'years' in cutoutparams and\
									'months' in cutoutparams and\
									all(os.path.isfile(self.datasetfn([y, m])) for y in range(cutoutparams['years'].start, cutoutparams['years'].stop+1) for m in range(cutoutparams['months'].start, cutoutparams['months'].stop+1) ):
				# All files are accounted for. Checking basic data and coverage

				if 'module' not in meta.attrs:
					raise TypeError('No module given in meta file of cutout.')
				# load dataset module based on file metadata

				self.dataset_module = sys.modules['geodata.datasets.' + meta.attrs['module']]
				cutoutparams['module'] = meta.attrs['module']

				logger.info("All cutout (%s, %s) files available.", name, cutout_dir)

				if {"xs", "ys"}.intersection(cutoutparams):
					# Passed some subsetting bounds
					self.meta = meta = self.get_meta_view(**cutoutparams)
					if not meta is None:
						# Subset is available
						self.prepared = True
						logger.info("Cutout subset prepared: %s", self)
					else:
						logger.info("Cutout subset not available: %s", self)
				else:
					# No subsetting of bounds. Keep full cutout
					self.prepared = True
					logger.info("Cutout prepared: %s", self)

			else:
				#   Not all files accounted for
				self.prepared = False
				logger.info("Cutout (%s, %s) not complete.", name, cutout_dir)

		if not self.prepared:
			# Still need to prepare cutout

			if 'module' not in cutoutparams:
				raise TypeError('Module is required to create cutout.')
			# load module from geodata library
			self.dataset_module = sys.modules['geodata.datasets.' + cutoutparams['module']]

			logger.info("Cutout (%s, %s) not found or incomplete.", name, cutout_dir)

			if {"xs", "ys", "years"}.difference(cutoutparams):
				raise TypeError("Arguments `xs`, `ys` and `years` need to be specified for a cutout.")

			if not meta is None:
				# if meta.nc exists, close and delete it
				meta.close()
				os.remove(self.datasetfn())

			## Main preparation call for metadata
			#    preparation.cutout_get_meta
			#    cutout.meta_data_config
			#    dataset_module.meta_data_config (e.g. prepare_meta_era5)
			self.meta = self.get_meta(**cutoutparams)

			# Ensure cutout directory exists
			if not os.path.isdir(self.cutout_dir):
				os.mkdir(self.cutout_dir)

			# Write meta file
			(self.meta_clean
				.unstack('year-month')
				.to_netcdf(self.datasetfn()))

	def datasetfn(self, *args):
		#    Link to dataset (default to meta.nc)
		#
		dataset = None

		if len(args) == 2:
			dataset = args
		elif len(args) == 1:
			dataset = args[0]
		else:
			dataset = None
		return os.path.join(self.cutout_dir, ("meta.nc"
											  if dataset is None
											  else "{}{:0>2}.nc".format(*dataset)))

	@property
	def meta_data_config(self):
		return dict(
			tasks_func=self.dataset_module.weather_data_config[self.config]['tasks_func'],
			prepare_func=self.dataset_module.weather_data_config[self.config]['meta_prepare_func'],
			template=self.dataset_module.weather_data_config[self.config]['template'],
			file_granularity=self.dataset_module.weather_data_config[self.config]['file_granularity']
		)
		#return self.dataset_module.meta_data_config
		## Step 2 - Change this to pull from
		## dict(
		##	prepare_func=self.weatherconfig['prepare_func'],
		##	template=self.weatherconfig['template']


	@property
	def weather_data_config(self):
		return self.dataset_module.weather_data_config

	@property
	def projection(self):
		return self.dataset_module.projection

	@property
	def coords(self):
		return self.meta.coords

	@property
	def meta_clean(self):
		# produce clean version of meta file for export to NetCDF (cannot export slices)
		meta = self.meta
		if meta.attrs.get('view', {}):
			view = {}
			for name, value in iteritems(meta.attrs.get('view', {})):
				view.update({name: [value.start, value.stop]})
			meta.attrs['view'] = str(view)
		return meta

	@property
	def shape(self):
		return len(self.coords["y"]), len(self.coords["x"])

	@property
	def extent(self):
		return (list(self.coords["x"].values[[0, -1]]) +
				list(self.coords["y"].values[[-1, 0]]))


	def grid_coordinates(self):
		xs, ys = np.meshgrid(self.coords["x"], self.coords["y"])
		return np.asarray((np.ravel(xs), np.ravel(ys))).T

	def grid_cells(self):
		from shapely.geometry import box
		coords = self.grid_coordinates()
		span = (coords[self.shape[1]+1] - coords[0]) / 2
		return [box(*c) for c in np.hstack((coords - span, coords + span))]

	def __repr__(self):
		yearmonths = self.coords['year-month'].to_index()
		return ('<Cutout {} x={:.2f}-{:.2f} y={:.2f}-{:.2f} time={}/{}-{}/{} {}prepared>'
				.format(self.name,
						self.coords['x'].values[0], self.coords['x'].values[-1],
						self.coords['y'].values[0], self.coords['y'].values[-1],
						yearmonths[0][0],  yearmonths[0][1],
						yearmonths[-1][0], yearmonths[-1][1],
						"" if self.prepared else "UN"))

	## Masking

	def add_mask(self, name, merged_mask = True, shape_mask = True):
		"""
		Add mask attribute to the cutout, from a previously saved mask objects.
		The masks will be coarsened to the same dimension with the cutout metadata in Xarray.

		name (str): name of the previously saved mask
		merged_mask (bool): if true, the program will try to 
			include the merged_mask from the mask object. Defaults to True.
		shape_mask (bool): if true, the program will try to 
			include the extracted dictionary of shape_mask from the mask object. Defaults to True.

		"""
		#make sure data is in correct format for coarsening
		xr_ds = self.meta.reset_coords(['lon', 'lat'], drop = True).rename({'x': 'lon', 'y': 'lat'})
		xr_ds = xr_ds.sortby(['lat', 'lon'])
		mask = load_mask(name, mask_dir = config.mask_dir)

		if not mask.merged_mask and not mask.shape_mask:
			raise ValueError(f"No mask found in {mask.name}. Please create a proper mask object first.")

		if mask.merged_mask and merged_mask:
			self.merged_mask = coarsen(mask.load_merged_xr(), xr_ds)
			logger.info("Cutout.merged_mask added.")

		if mask.shape_mask and shape_mask:
			self.shape_mask = {k: coarsen(v, xr_ds) for k, v in mask.load_shape_xr().items()}
			logger.info("Cutout.shape_mask added.")

	def add_grid_area(self, axis = ("lat", "lon"), adjust_coord = True):
		"""
		Add attribute 'area' to the cutout containing area for each grid cell 
		in the cutout metedata Xarray.

		axis (tuple): axis to include in the result xarray dataset. 
			Defaults to ("lat", "lon").
		adjust_coord (bool): sort the data by latitude and longitude values if true. 
			Defaults to True.

		"""
		xr_ds = self.meta.reset_coords(['lon', 'lat'], drop = True).rename({'x': 'lon', 'y': 'lat'})
		area_arr = np.zeros((xr_ds.lat.shape[0], xr_ds.lon.shape[0]))
		lat_diff = np.abs((xr_ds.lat[1].values - xr_ds.lat[0].values))
		for i, lat in enumerate(xr_ds.lat.values):
			lat_bottom = lat - lat_diff/ 2
			lat_top = lat + lat_diff/ 2
			#calculate the area for grid cells with same latitude
			area_arr[i] = np.round(calc_grid_area([
							(xr_ds.lon.values[0], lat_top), 
							(xr_ds.lon.values[0], lat_bottom),
							(xr_ds.lon.values[1], lat_bottom), 
							(xr_ds.lon.values[1], lat_top)]), 2)
		if axis == ("lat", "lon"):
			xr_ds = xr_ds.assign({"area": (axis, area_arr)})
		elif axis == ("time", "lat", "lon"):
			xr_ds = xr_ds.assign({"area": (axis, [area_arr])})

		if adjust_coord:
			if (('lon' not in xr_ds.coords) and ('lat' not in xr_ds.coords)
				and ('x' in xr_ds.coords and 'y' in xr_ds.coords)):
				xr_ds = xr_ds.rename({'x': 'lon', 'y': 'lat'})
			xr_ds = xr_ds.sortby(['lat', 'lon'])

		self.area = xr_ds
	
	def mask(self, dataset, true_area = True, 
					merged_mask = True, shape_mask = True):
		"""
		Mask a converted Xarray dataSet from cutout with previously added mask attribute
		with optional area inclusion, and return a dictionary of xarray Dataset.

		The program will search for 'merged_mask' and 'shape_mask' attributes in the 
		cutout object, these Xarray data can be generate through 'add_mask', unless the user 
		specify 'merged_mask = False' or 'shape_mask = False', the masks in shape_mask 
		will have the same key in the dictionary returned, and the mask for merged_mask will 
		have the key name "merged_mask".

		dataset (Xarray.DataSet): 
		true_area (bool): if the returned masks will have the area variable. Defaults to True.
		merged_mask (bool): if true, the program will try to 
			include the merged_mask from the cutout object. Defaults to True.
		shape_mask (bool): if true, the program will try to 
			include the extracted dictionary of shape_mask from the cutout object. Defaults to True.

		Returns: (dict): a dictionary with name keys and xarray DataSet with mask variable as values.
		
		"""
		axis = ("lat", "lon")

		if self.merged_mask is None and self.shape_mask is None:
			raise ValueError(f"No mask found in cutout. Please add masks with self.add_mask()")

		if self.area is None and true_area == True:
			raise ValueError("No area data found. Please call self.add_grid_area() or set true_area to False.")

		res = {}
		
		if self.merged_mask is not None and merged_mask:
			ds = dataset.assign({"mask": (axis, self.merged_mask.data[0])})
			if true_area:
				ds = ds.assign({"area": (axis, self.area['area'].data)})
			res['merged_mask'] = ds
			logger.info("merged_mask combined with dataset. ")

		if self.shape_mask is not None and shape_mask:
			for key, val in self.shape_mask.items():
				ds = dataset.assign({"mask": (axis, val.data[0])})
				if true_area:
					ds = ds.assign({"area": (axis, self.area['area'])})
				res[key] = ds
		logger.info("shape_mask combined with dataset. ")

		return res
		
	

	## Preparation functions

	get_meta = cutout_get_meta    # preparation.cutout_get_meta

	get_meta_view = cutout_get_meta_view  # preparation.cutout_get_meta_view

	prepare = cutout_prepare      # preparation.cutout_prepare

	produce_specific_dataseries = cutout_produce_specific_dataseries

	## Conversion and aggregation functions

	convert_cutout = convert_cutout

	heat_demand = heat_demand

	temperature = temperature

	soil_temperature = soil_temperature

	solar_thermal = solar_thermal

	wind = wind

	pv = pv



def _find_intercept(list1, list2, start, threshold = 0):
	'''
	find_intercept is a helper function to find the best start point for doing coarsening 
	in order to make the coordinates of the coarsen as close to the target as possible.
	'''
	min_pos = 0
	min_res = 0
	init = 0
	for i in range(len(list1)-start):
		resid = ((list1[start+i] -list2[0]) % (list2[1] - list2[0])).values.tolist()
		if i == 0:
			init = resid
		if resid <= threshold:
			return i
		if resid > min_res:
			min_res = resid
			min_pos = i
		else:
			min_res = resid
			break
	if min_res == init:
		return 0
	else:
		return i #type: ignore

def coarsen(ori, tar, threshold = 0, func = 'mean'):
	'''
	This function will reindex the original xarray dataset according to the coordiantes of the target.
	
	There might be a bias for lattitudes and longitudes. The bias are normally within 0.01 degrees. 
	In order to not lose too much data, a threshold for bias in degree could be given. 
	When threshold = 0, it means that the function is going to find the best place with smallest bias.
	
	'''
	lat_multiple = round(((tar.lat[1] - tar.lat[0]) / (ori.lat[1] - ori.lat[0])).values.tolist())
	lon_multiple = round(((tar.lon[1] - tar.lon[0]) / (ori.lon[1] - ori.lon[0])).values.tolist())
	
	lat_start = _find_intercept(ori.lat, tar.lat, (lat_multiple - 1) //  2) 
	lon_start = _find_intercept(ori.lon, tar.lon, (lon_multiple - 1) //  2)
	
	if func == 'mean':
		coarsen = ori.isel(lat = slice(lat_start, None), lon = slice(lon_start, None)).coarsen(
			dim = {'lat': lat_multiple, 'lon' : lon_multiple}, 
			side = {'lat': 'left', 'lon': 'left'}, 
			boundary = 'pad', ).mean()
	elif func == 'sum':
		coarsen = ori.isel(lat = slice(lat_start, None), lon = slice(lon_start, None)).coarsen(
			dim = {'lat': lat_multiple, 'lon' : lon_multiple}, 
			side = {'lat': 'left', 'lon': 'left'}, 
			boundary = 'pad', ).sum()
	out = coarsen.reindex_like(tar, method = 'nearest')
	return out	

def calc_grid_area(lis_lats_lons):
	"""
	Calculate area in km^2 for a grid cell given lats and lon border, with help from:
	https://stackoverflow.com/questions/4681737/how-to-calculate-the-area-of-a-polygon-on-the-earths-surface-using-python

	"""
	lons, lats = zip(*lis_lats_lons)
	ll = list(set(lats))[::-1]
	var = []
	for i in range(len(ll)):
		var.append('lat_' + str(i+1))
	st = ""
	for v, l in zip(var,ll):
		st = st + str(v) + "=" + str(l) +" "+ "+"
	st = st +"lat_0="+ str(np.mean(ll)) + " "+ "+" + "lon_0" +"=" + str(np.mean(lons))
	tx = "+proj=aea +" + st
	pa = pyproj.Proj(tx)

	x, y = pa(lons, lats)
	cop = {"type": "Polygon", "coordinates": [zip(x, y)]}

	return shapely.geometry.shape(cop).area/1000000

def calc_shp_area(shp, shp_projection = '+proj=latlon'):
	"""calculate area in km^2 of the shapes for each shp object"""
	temp_shape =  shapely.ops.transform(partial(
            pyproj.transform,
            pyproj.Proj(shp_projection),
            pyproj.Proj(
                proj='aea',
                lat_1 = shp.bounds[1],
                lat_2 = shp.bounds[3]
            )
        ),
    shp)
	return temp_shape.area / 1000000
