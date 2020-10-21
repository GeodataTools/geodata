## Copyright 2020 Michael Davidson (UCSD), William Honaker. Copyright 2016-2017 Gorm Andresen (Aarhus University), Jonas Hoersch (FIAS), Tom Brown (FIAS)

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
from six import string_types, itervalues, iteritems

import logging
logger = logging.getLogger(__name__)

from . import config, datasets

from .convert import (convert_and_aggregate, heat_demand, hydro, temperature,
					  wind, pv, runoff, solar_thermal, soil_temperature)
from .preparation import (cutout_do_task, cutout_prepare,
						  cutout_produce_specific_dataseries,
						  cutout_get_meta, cutout_get_meta_view)
from .gis import compute_indicatormatrix

class Cutout(object):
	def __init__(self, name=None, cutout_dir=config.cutout_dir, **cutoutparams):
		self.name = name
		self.cutout_dir = os.path.join(cutout_dir, name)
		self.prepared = False
		self.meta_append = 0
		self.config = cutoutparams.pop('weather_data_config')
		self.meta = meta = None


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

	def indicatormatrix(self, shapes, shapes_proj='latlong'):
		return compute_indicatormatrix(self.grid_cells(), shapes, self.projection, shapes_proj)

	## Preparation functions

	get_meta = cutout_get_meta    # preparation.cutout_get_meta

	get_meta_view = cutout_get_meta_view  # preparation.cutout_get_meta_view

	prepare = cutout_prepare      # preparation.cutout_prepare

	produce_specific_dataseries = cutout_produce_specific_dataseries

	## Conversion and aggregation functions

	convert_and_aggregate = convert_and_aggregate

	heat_demand = heat_demand

	temperature = temperature

	soil_temperature = soil_temperature

	solar_thermal = solar_thermal

	wind = wind

	pv = pv

	runoff = runoff

	hydro = hydro
