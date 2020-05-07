## dataset.py
#	dataset class added to Atlite (https://github.com/PyPSA/atlite/) by Michael Davidson (UCSD)
#	 - downloading and verifying local datasets downloaded (in config.py)

"""
Renewable Energy Atlas Lite (Atlite)

Light-weight version of Aarhus RE Atlas for converting weather data to power systems data
"""

from __future__ import absolute_import

import xarray as xr
import numpy as np
import os, sys, glob, requests, shutil
from requests.exceptions import HTTPError
from tempfile import mkstemp
from six import string_types, itervalues, iteritems
from calendar import monthrange

import logging
logger = logging.getLogger(__name__)

from . import config, datasets


class Dataset(object):
	def __init__(self, **datasetparams):
		if 'module' not in datasetparams:
			raise ValueError("`module` needs to be specified")
		if 'weather_data_config' not in datasetparams:
			raise ValueError("`config` needs to be specified")
		self.module = datasetparams.pop('module')
		self.config = datasetparams.pop('weather_data_config')
		# load module from geodata library
		self.dataset_module = sys.modules['geodata.datasets.' + self.module]
		
		self.weatherconfig = self.weather_data_config[self.config]

		if 'datadir' in datasetparams:
			logger.info("Manual data directory entry not supported. Change in config.py file.")
		# 	self.datadir = datasetparams.pop('datadir')
		# else:
		self.datadir = self.dataset_module.datadir

		if not "years" in datasetparams:
			raise ValueError("`years` need to be specified")
		self.years = years = datasetparams['years']

		if not "months" in datasetparams:
			logger.info("No months specified, defaulting to 1-12")
			self.months = months = slice(1,12)
		else:
			self.months = months = datasetparams['months']

		self.prepared = False
		self.toDownload = []
		self.savedFiles = []
		self.downloadedFiles = []
		self.totalFiles = []
		incomplete_count = 0

		if 'bounds' in datasetparams:
			x1, y1, x2, y2 = datasetparams.pop('bounds')
			datasetparams.update(xs=slice(x1, x2),
								ys=slice(y2, y1))

		if os.path.isdir(self.datadir):
			# Directory for dataset exists
			logger.info("Directory %s found, checking for completeness.", self.datadir)
			self.prepared = True
			check_complete = True
		else:
			logger.info("Directory %s not found.", self.datadir)
			check_complete = False

		step = years.step if years.step else 1
		yrs = range(years.start, years.stop+step, step)
		step = months.step if months.step else 1
		mos = range(months.start, months.stop+step, step)

		if self.weatherconfig['time_period'] == 'daily':
			# Separate files for each day (eg MERRA)
			# Check for complete set of files of year, month, day
			mo_tuples = [(yr,mo,monthrange(yr,mo)[1]) for yr in yrs for mo in mos]

			for mo_tuple in mo_tuples:
				# format: (yr, mo, number_days_in_month)
				yr, mo, nodays = mo_tuple
				for day in range(1, nodays+1, 1):
					filename = self.datasetfn(self.weatherconfig['fn'], yr, mo, day)
					self.totalFiles.append((self.config, filename))
					if not os.path.isfile( filename ):
						self.prepared = False
						if check_complete:
							logger.info("File `%s` not found!", filename)
							incomplete_count += 1
						self.toDownload.append((self.config, filename, self.datasetfn(self.weatherconfig['url'], yr, mo, day)))
					else:
						self.downloadedFiles.append((self.config, filename))

		elif self.weatherconfig['time_period'] == 'monthly':
			# Monthly files (eg ERA5)
			mo_tuples = [(yr,mo) for yr in yrs for mo in mos]
			for mo_tuple in mo_tuples:
				yr, mo = mo_tuple
				filename = self.datasetfn(self.weatherconfig['fn'], yr, mo)
				self.totalFiles.append((self.config, filename))
				if not os.path.isfile( filename ):
					self.prepared = False
					if check_complete:
						logger.info("File `%s` not found!", filename)
						incomplete_count += 1
					self.toDownload.append((self.config, filename, self.datasetfn(self.weatherconfig['url'], yr, mo)))
				else:
					self.downloadedFiles.append((self.config, filename))

		 # removed unneeded check for one file section

		if not self.prepared:

			if {"xs", "ys"}.difference(datasetparams):
				logger.warn("Arguments `xs` and `ys` not used in preparing dataset. Defaulting to global.")

			logger.info(f'{incomplete_count} files not completed.')
			## Main preparation call for metadata
			#	preparation.cutout_get_meta
			#	cutout.meta_data_config
			#	dataset_module.meta_data_config (e.g. prepare_meta_era5)
			# self.meta = self.get_meta(**datasetparams)
			return None
		else:
			logger.info("Directory complete.")
			return None

	def datasetfn(self, fn, *args):
		# construct file name from fn template (cf weather_data_config) and args (yr, mo, day)
		if len(args) == 3:
			dataset = dict(year=args[0], month=args[1], day=args[2])
		elif len(args) == 2:
			dataset = dict(year=args[0], month=args[1])
		else:
			return False
		if self.dataset_module.spinup_var:
			spinup = self.dataset_module.spinup_year(dataset['year'])
			dataset.update({'spinup': spinup})
			
		return fn.format_map(dataset)

	def get_data(self, trim=False, testing=False, wind=True, solar=True):
		"""Download data routine
		# By default, keep variables related to wind (True) and solar (True)
		#
		#	Parameters
		#	---------
		#	trim: boolean
		#		Run trim_variables function following each download
		"""
		# Changes 5/5/2020
		# Testing parameter removed in favor of future functionality around daily downloads.
		# Get data now happens on a per config basis, as dataset is also defined on a per-config basis.
		# Dataset now updates to "prepared" status following completion of download, without having to
		# redfine object.

		self.savedFiles = []

		count = 0
		for f in self.toDownload:
			print(f)

			# Make the directory if not exists:
			os.makedirs(os.path.dirname(f[1]), exist_ok=True)
			result = requests.get(f[2])
			try:
					result.raise_for_status()
					fout = open(f[1],'wb')
					fout.write(result.content)
					fout.close()
					self.savedFiles.append((f[0], f[1])) # What is saved files being used for?
					if trim:
						self.trim_variables( fn = [(f[0], f[1])], wind = wind, solar=solar )
			except HTTPError as http_err:
					logger.warn(f'HTTP error occurred: {http_err}')  # Python 3.6
			except Exception as err:
					logger.warn(f'Other error occurred: {err}')  # Python 3.6
					# logger.warn('requests.get() returned an error code '+str(result.status_code))
			count += 1

		if self.savedFiles == self.totalFiles: 
			self.prepared = True	
		# Added functionality to auto update dataset object
		# to prepared = True upon completion of download.



#		weather_data = []
#		if wind:
#			weather_data.extend(self.dataset_module.wind_files)
#		if solar:
#			weather_data.extend(self.dataset_module.solar_files)
#		weather_data = list(set(weather_data))
#
#		# Loop through files identified as missing in constructor
#		count = 0
#		for f in self.toDownload:
#			if testing and (count == 1):
#				# (for testing purposes) download only the first file
#				continue
#			if f[0] in weather_data:
#				# File is in list of datasets we want to download
#				print(f)
#
#				# Make the directory if not exists:
#				os.makedirs(os.path.dirname(f[1]), exist_ok=True)

#				result = requests.get(f[2])
#				try:
#					result.raise_for_status()
#					fout = open(f[1],'wb')
#					fout.write(result.content)
#					fout.close()
#					self.savedFiles.append((f[0], f[1]))
#					if trim:
#						self.trim_variables( fn = [(f[0], f[1])], wind = wind, solar=solar )
#				except HTTPError as http_err:
#					logger.warn(f'HTTP error occurred: {http_err}')  # Python 3.6
#				except Exception as err:
#					logger.warn(f'Other error occurred: {err}')  # Python 3.6
					# logger.warn('requests.get() returned an error code '+str(result.status_code))
#			count += 1 */

	def trim_variables(self, fn = None, downloadedfiles = False, wind=True, solar=True):
		""" Reduce size of file by trimming variables in file
		# 	By default, keep variables related to wind (True) and solar (True)
		#
		#	Parameters
		#	---------
		#	fn: array of tuples
		#		if present: use this array of files
		#	downloadedfiles: boolean
		#		if True: use downloadedFiles array (of files already on system during dataset construction)
		#		if False: use savedFiles array (of files downloaded from get_data())
		#
		# TODO: create options for non NetCDF4 files (eg pynio)
		"""

		if downloadedfiles:
			filelist = self.downloadedFiles
		elif not fn is None:
			filelist = fn
		else:
			filelist = self.savedFiles

		for d, f in filelist:
			# d = filetype (in weather_data_config), f = filename

			# Construct list of variables to keep
			vars = []
			if wind and d in self.dataset_module.wind_files:
				vars.extend(self.weather_data_config[d]['variables'])
			if solar and d in self.dataset_module.solar_files:
				vars.extend(self.weather_data_config[d]['variables'])
			vars = list(set(vars))

			with xr.open_dataset(f) as ds:
				# vars to lower case
				var_rename = dict((v, v.lower()) for v in list(ds.data_vars))
				ds = ds.rename(var_rename)
				ds = ds[vars]

				# Save to temp file
				fd, target = mkstemp(suffix='.nc4')
				os.close(fd)
				ds.to_netcdf(target)

			#
			# newf = os.path.join(os.path.dirname(f),'trim/',os.path.split(f)[1])
			# # Make the directory if not exists:
			# os.makedirs(os.path.dirname(newf), exist_ok=True)
			# ds.to_netcdf(newf)

			# Move temp file
			shutil.move(target,f)



	def set_saved_files(self):
		self.savedFiles = [('surface_flux', '/Users/michd/Documents/GEODATA/data/merra2/2011/01/MERRA2_400.tavg1_2d_flx_Nx.20110101.nc4')]


	@property
	def meta_data_config(self):
		return self.dataset_module.meta_data_config

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
		return ('<Dataset {} years={}-{} months={}-{} datadir={} {}>'
				.format(self.module,
						self.years.start, self.years.stop,
						self.months.start, self.months.stop,
						self.datadir,
						"Prepared" if self.prepared else "Unprepared"))

	def indicatormatrix(self, shapes, shapes_proj='latlong'):
		return compute_indicatormatrix(self.grid_cells(), shapes, self.projection, shapes_proj)