## Copyright 2016-2017 Jonas Hoersch (FIAS), Tom Brown (FIAS), Markus Schlott
## (FIAS)

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

import os
import glob
import pandas as pd
import numpy as np
import xarray as xr
import shutil
from six.moves import range
from contextlib import contextmanager
from tempfile import mkstemp
import logging
logger = logging.getLogger(__name__)

from ..config import era5_dir

datadir = era5_dir

try:
	import cdsapi
	has_cdsapi = True
except ImportError:
	has_cdsapi = False

# Model and Projection Settings
projection = 'latlong'

# contextmanager (contextlib): ensures that at end of all with statements, the resource is closed
@contextmanager
def _get_data(target=None, product='reanalysis-era5-single-levels', chunks=None, **updates):
	"""
	Check local folders for ERA5
	If necessary: download ERA5 data from the Climate Data Store (CDS)
	"""

	## Check if local copy
	# filename = first characters of each variable
	f = "".join([v[0] for v in updates['variable']])
	if target is None:
		target = era5_dir

	fn=os.path.join(era5_dir, '{year}/{month:0>2}/{f}.nc'.format(year=updates['year'], month=updates['month'],f=f))

	if os.path.isfile(fn):
	#	Local file exists
		with xr.open_dataset(fn, chunks=chunks) as ds:
			if {'area'}.issubset(updates):
			# find subset of area
				yf, x0, y0, xf = updates['area']			# North, West, South, East
				ds = ds.where((ds.latitude >= y0) & (ds.latitude <= yf) & (ds.longitude >= x0) & (ds.longitude <= xf), drop=True)

			yield ds
	else:
	#	Download new file

		# Make the directory if not exists
		os.makedirs(os.path.dirname(fn), exist_ok=True)

		if not has_cdsapi:
			raise RuntimeError(
				"Need installed cdsapi python package available from "
				"https://cds.climate.copernicus.eu/api-how-to"
			)

		# Default request
		request = {
			'product_type':'reanalysis',
			'format':'netcdf',
			'day':[
				'01','02','03','04','05','06','07','08','09','10','11','12',
				'13','14','15','16','17','18','19','20','21','22','23','24',
				'25','26','27','28','29','30','31'
			],
			'time':[
				'00:00','01:00','02:00','03:00','04:00','05:00',
				'06:00','07:00','08:00','09:00','10:00','11:00',
				'12:00','13:00','14:00','15:00','16:00','17:00',
				'18:00','19:00','20:00','21:00','22:00','23:00'
			]
		}
		request.update(updates)

		assert {'year', 'month', 'variable'}.issubset(request), "Need to specify at least 'variable', 'year' and 'month'"

		# main request to API
		result = cdsapi.Client().retrieve(
			product,
			request
		)

		# if target is None:
		# 	# make a temp file
		#     fd, target = mkstemp(suffix='.nc')
		#     os.close(fd)
		logger.info("Downloading request for {} variables to {}".format(len(request['variable']), target))
		result.download(fn)

		with xr.open_dataset(fn, chunks=chunks) as ds:
			yield ds

		# os.unlink(target)

def _add_height(ds):
	"""Convert geopotential 'z' to geopotential height following [1]

	References
	----------
	[1] ERA5: surface elevation and orography, retrieved: 10.02.2019
	https://confluence.ecmwf.int/display/CKB/ERA5%3A+surface+elevation+and+orography

	"""
	g0 = 9.80665
	z = ds['z']
	if 'time' in z.coords:
		z = z.isel(time=0, drop=True)
	ds['height'] = z/g0
	ds = ds.drop('z')
	return ds

def _area(xs, ys):
	# Return array with bounding coordinates
	# North, West, South, East.
	return [ys.start, xs.start, ys.stop, xs.stop]

def _rename_and_clean_coords(ds, add_lon_lat=True):
	"""Rename 'longitude' and 'latitude' columns to 'x' and 'y'

	Optionally (add_lon_lat, default:True) preserves latitude and longitude columns as 'lat' and 'lon'.
	"""

	ds = ds.rename({'lon': 'x', 'lat': 'y'})
	if add_lon_lat:
		ds = ds.assign_coords(lon=ds.coords['x'], lat=ds.coords['y'])
	return ds

def convert_and_subset_lons_lats_era5(ds, xs, ys):
	# Rename geographic dimensions to x,y
	# Subset x,y according to xs, ys

	if not isinstance(xs, slice):
		first, second, last = np.asarray(xs)[[0,1,-1]]
		xs = slice(first - 0.1*(second - first), last + 0.1*(second - first))
	if not isinstance(ys, slice):
		first, second, last = np.asarray(ys)[[0,1,-1]]
		ys = slice(first - 0.1*(second - first), last + 0.1*(second - first))

	ds = ds.sel(lat=ys)

	# Lons should go from -180. to +180.
	if len(ds.coords['lon'].sel(lon=slice(xs.start + 360., xs.stop + 360.))):
		ds = xr.concat([ds.sel(lon=slice(xs.start + 360., xs.stop + 360.)),
						ds.sel(lon=xs)],
					   dim="lon")
		ds = ds.assign_coords(lon=np.where(ds.coords['lon'].values <= 180,
											 ds.coords['lon'].values,
											 ds.coords['lon'].values - 360.))
	else:
		ds = ds.sel(lon=xs)

	ds = ds.rename({'lon': 'x', 'lat': 'y'})
	ds = ds.assign_coords(lon=ds.coords['x'], lat=ds.coords['y'])
	return ds


def subset_x_y_era5(ds, xs, ys):
	# Subset x,y according to xs, ys

	if not isinstance(xs, slice):
		first, second, last = np.asarray(xs)[[0,1,-1]]
		xs = slice(first - 0.1*(second - first), last + 0.1*(second - first))
	if not isinstance(ys, slice):
		first, second, last = np.asarray(ys)[[0,1,-1]]
		ys = slice(first - 0.1*(second - first), last + 0.1*(second - first))

	ds = ds.sel(y=ys)
	ds = ds.sel(x=xs)

	return ds

def prepare_meta_era5(xs, ys, year, month, template, module, **kwargs):
	# Reference of the quantities
	# https://confluence.ecmwf.int/display/CKB/ERA5+data+documentation
	# Geopotential is aka Orography in the CDS:
	# https://confluence.ecmwf.int/pages/viewpage.action?pageId=78296105

	fns = glob.iglob(template.format(year=year, month=month))
	with xr.open_mfdataset(fns, combine='by_coords') as ds:
		ds = ds.coords.to_dataset()
		ds = convert_and_subset_lons_lats_era5(ds, xs, ys)
		meta = ds.load()

	return meta


def prepare_for_sarah(year, month, xs, ys, dx, dy, chunks=None):
	area = _area(xs, ys)
	grid = [dx, dy]

	with _get_data(area=area, grid=grid, year=year, month=month,
				   variable=['2m_temperature',
							 'toa_incident_solar_radiation',
							 'surface_solar_radiation_downwards',
							 'surface_net_solar_radiation'],
				   chunks=chunks) as ds:
		ds = _rename_and_clean_coords(ds, add_lon_lat=False)

		ds = ds.rename({'t2m': 'temperature'})
		ds = ds.rename({'tisr': 'influx_toa'})

		logger.debug("Calculating albedo")
		ds['albedo'] = (((ds['ssrd'] - ds['ssr'])/ds['ssrd']).fillna(0.)
						.assign_attrs(units='(0 - 1)', long_name='Albedo'))
		ds = ds.drop(['ssrd', 'ssr'])

		logger.debug("Fixing units of influx_toa")
		# Convert from energy to power J m**-2 -> W m**-2
		ds['influx_toa'] /= 60.*60.
		ds['influx_toa'].attrs['units'] = 'W m**-2'

		logger.debug("Yielding ERA5 results")
		yield ds.chunk(chunks)
		logger.debug("Cleaning up ERA5")

def prepare_month_era5(fn, year, month, xs, ys):

	# Reference of the quantities
	# https://confluence.ecmwf.int/display/CKB/ERA5+data+documentation
	# (shortName) | (name)                                      | (paramId)
	# tisr        | TOA incident solar radiation                | 212
	# ssrd        | Surface Solar Rad Downwards                 | 169
	# ssr         | Surface net Solar Radiation                 | 176
	# fdir        | Total sky direct solar radiation at surface | 228021
	# ro          | Runoff                                      | 205
	# 2t          | 2 metre temperature                         | 167
	# sp          | Surface pressure                            | 134
	# stl4        | Soil temperature level 4                    | 236
	# fsr         | Forecast surface roughnes                   | 244

	if not os.path.isfile(fn):
		return None
	with xr.open_dataset(fn) as ds:
		logger.info(f'Opening `{fn}`')
		ds = _rename_and_clean_coords(ds)
		ds = subset_x_y_era5(ds, xs, ys)
		yield (year, month), ds


def tasks_monthly_era5(xs, ys, yearmonths, prepare_func, **meta_attrs):
	if not isinstance(xs, slice):
		xs = slice(*xs.values[[0, -1]])
	if not isinstance(ys, slice):
		ys = slice(*ys.values[[0, -1]])
	fn = meta_attrs['fn']

	logger.info(yearmonths)
	logger.info([(year, month) for year, month in yearmonths])

	return [
		dict(
			prepare_func=prepare_func, 
			xs=xs, 
			ys=ys, 
			year=year, 
			month=month,
			fn=fn.format(year=year, month=month)
		)
		for year, month in yearmonths
	]

weather_data_config = {
	'era5_monthly': dict(
		file_granularity="monthly",
		tasks_func=tasks_monthly_era5,
		meta_prepare_func=prepare_meta_era5,
		prepare_func=prepare_month_era5,
		template=os.path.join(era5_dir, '{year}/{month:0>2}/*.nc'),
		fn = os.path.join(era5_dir, '{year}/{month:0>2}/112rsssstt.nc'),
		product='reanalysis-era5-single-levels',
		variables=[
					   '100m_u_component_of_wind',
					   '100m_v_component_of_wind',
					   '2m_temperature',
					   'runoff',
					   'soil_temperature_level_4',
					   'surface_net_solar_radiation',
					   'surface_pressure',
					   'surface_solar_radiation_downwards',
					   'toa_incident_solar_radiation',
					   'total_sky_direct_solar_radiation_at_surface'
				   ]
		)
}

#meta_data_config = dict(prepare_func=prepare_meta_era5)

# No separate files for each day (would be coded in weather_data_config list, see merra2.py)
daily_files = False

# Latitude stored south to north (ie forward, = True) or north to south
lat_direction = False

# Spinup variable (necessary for MERRA)
spinup_var = False
