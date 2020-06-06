## merra2.py
#	MERRA-2 dataset routine added to Atlite (https://github.com/PyPSA/atlite/) by Michael Davidson (UCSD)

#

"""
Renewable Energy Atlas Lite (Atlite)

Light-weight version of Aarhus RE Atlas for converting weather data to power systems data
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
from calendar import monthrange

import logging
logger = logging.getLogger(__name__)

from ..config import merra2_dir

datadir = merra2_dir

# Model and Projection Settings
projection = 'latlong'

def convert_and_subset_lons_lats_merra2(ds, xs, ys):
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

def subset_x_y_merra2(ds, xs, ys):
	# Subset x,y according to xs, ys
	# Assumes convert_and_subset_lons_lats_merra2 already run

	if not isinstance(xs, slice):
		first, second, last = np.asarray(xs)[[0,1,-1]]
		xs = slice(first - 0.1*(second - first), last + 0.1*(second - first))
	if not isinstance(ys, slice):
		first, second, last = np.asarray(ys)[[0,1,-1]]
		ys = slice(first - 0.1*(second - first), last + 0.1*(second - first))

	ds = ds.sel(y=ys)
	ds = ds.sel(x=xs)

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


def prepare_meta_merra2(xs, ys, year, month, template, module, **params):
	#	Load dataset into metadata

	# fn = next(glob.iglob(template.format(year=year, month=month)))
	# with xr.open_dataset(fn) as ds:
	# 	ds = ds.coords.to_dataset()
	# 	ds = convert_and_subset_lons_lats_merra2(ds, xs, ys)
	# 	meta = ds.load()

	# Set spinup variable (see MERRA2 documentation, p. 13)
	spinup = spinup_year(year)

	fns = glob.iglob(template.format(year=year, month=month, spinup=spinup))
	with xr.open_mfdataset(fns, combine='by_coords') as ds:
		ds = ds.coords.to_dataset()
		ds = convert_and_subset_lons_lats_merra2(ds, xs, ys)
		meta = ds.load()

	# Any time adjustments?
			# t = pd.Timestamp(year=year, month=month, day=1)
			# ds['time'] = pd.date_range(t, t + pd.DateOffset(months=1),
			# 						   freq='1h', closed='left')

	return meta


def prepare_month_surface_flux(fn, year, month, xs, ys):
	if not os.path.isfile(fn):
		return None
	with xr.open_dataset(fn) as ds:
		logger.info(f'Opening `{fn}`')
		# logger.info("Cutout dims: %s", ds.dims)
		# logger.info("Cutout coords: %s", ds.coords)

		ds = _rename_and_clean_coords(ds)

		# logger.info("Cutout dims: %s", ds.dims)
		# logger.info("Cutout coords: %s", ds.coords)

		ds = subset_x_y_merra2(ds, xs, ys)

		# logger.info("Cutout dims: %s", ds.dims)
		# logger.info("Cutout coords: %s", ds.coords)


		# some variable renaming
		try:
		#	z0m=roughness
		#	wind variables not in wndXXm format
			ds = ds.rename({'z0m': 'roughness'})
		except Exception as e:
			logger.warn(f'Unable to rename variables in `{fn}`. Exception: {e}')

		ds['wndlml'] = (np.sqrt(ds['ulml']**2 + ds['vlml']**2)
						.assign_attrs(units=ds['ulml'].attrs['units'],
									long_name="LML wind speed"))
		if ds['tlml']:
			ds['temperature'] = ds['tlml']

		yield (year, month), ds

def prepare_dailymeans_surface_flux(fn, year, month, xs, ys):
	if not os.path.isfile(fn):
		return None
	with xr.open_dataset(fn) as ds:
		logger.info(f'Opening `{fn}`')
		# logger.info("Cutout dims: %s", ds.dims)
		# logger.info("Cutout coords: %s", ds.coords)

		ds = _rename_and_clean_coords(ds)

		# logger.info("Cutout dims: %s", ds.dims)
		# logger.info("Cutout coords: %s", ds.coords)

		ds = subset_x_y_merra2(ds, xs, ys)

		# logger.info("Cutout dims: %s", ds.dims)
		# logger.info("Cutout coords: %s", ds.coords)

				# some variable renaming
		try:
		#	z0m=roughness
		#	wind variables not in wndXXm format
			ds = ds.rename({
				'T2MMEAN': 'temperature',
				'TPRECMAX': 'precipitation'
				})
		except Exception as e:
			logger.warn(f'Unable to rename variables in `{fn}`. Exception: {e}')

		#['HOURNORAIN', 'T2MMAX', 'T2MMEAN', 'T2MMIN', 'TPRECMAX']

		yield (year, month), ds



## TODO def prepare_month_radiation
		# with np.errstate(divide='ignore', invalid='ignore'):
		# 	ds['albedo'] = (((ds['ssrd'] - ds['ssr'])/ds['ssrd']).fillna(0.)
		# 					.assign_attrs(units='(0 - 1)', long_name='Albedo'))
		# ds['influx_diffuse'] = ((ds['ssrd'] - ds['influx_direct'])
		# 						.assign_attrs(units='J m**-2',
		# 									long_name='Surface diffuse solar radiation downwards'))
		# ds = ds.drop(['ssrd', 'ssr'])
		#
		# # Convert from energy to power J m**-2 -> W m**-2 and clip negative fluxes
		# for a in ('influx_direct', 'influx_diffuse', 'influx_toa'):
		# 	ds[a] = ds[a].clip(min=0.) / (60.*60.)
		# 	ds[a].attrs['units'] = 'W m**-2'


def tasks_daily_merra2(xs, ys, yearmonths, prepare_func, **meta_attrs):
	if not isinstance(xs, slice):
		xs = slice(*xs.values[[0, -1]])
	if not isinstance(ys, slice):
		ys = slice(*ys.values[[0, -1]])
	fn = meta_attrs['fn']

	logger.info(yearmonths)
	logger.info([(year, month, day) for year, month in yearmonths for day in range(1, monthrange(year,month)[1]+1, 1)])

	return [dict(prepare_func=prepare_func,
				 xs=xs, ys=ys,
				 year=year, month=month,
				 fn=fn.format(year=year, month=month, day=day, spinup=spinup_year(year)) )
				 for year, month in yearmonths for day in range(1, monthrange(year,month)[1]+1, 1)]

def tasks_monthly_merra2(xs, ys, yearmonths, prepare_func, **meta_attrs):
	if not isinstance(xs, slice):
		xs = slice(*xs.values[[0, -1]])
	if not isinstance(ys, slice):
		ys = slice(*ys.values[[0, -1]])
	fn = meta_attrs['fn']

	logger.info(yearmonths)
	logger.info([(year, month) for year, month in yearmonths])

	return [dict(prepare_func=prepare_func,
				 xs=xs, ys=ys,
				 year=year, month=month,
				 fn=fn.format(year=year, month=month, spinup=spinup_year(year)) )
				 for year, month in yearmonths]



weather_data_config = {
#	Single file contains all wind variables (≠ ncep)
#	MERRA2 has additional label for spinup decade--eg 300, 400--that must be calculated via spinup_year(year) before downloading
# 	https://goldsmr4.gesdisc.eosdis.nasa.gov/data/MERRA2/M2T1NXFLX.5.12.4/2015/01/MERRA2_400.tavg1_2d_flx_Nx.20150101.nc4
	'surface_flux_hourly': dict(
		file_granularity="daily",
		tasks_func=tasks_daily_merra2,
		meta_prepare_func=prepare_meta_merra2,
		prepare_func=prepare_month_surface_flux,
		template=os.path.join(merra2_dir, '{year}/{month:0>2}/MERRA2_*.tavg1_2d_flx_Nx.*.nc4'),
		url = 'https://goldsmr4.gesdisc.eosdis.nasa.gov/data/MERRA2/M2T1NXFLX.5.12.4/{year}/{month:0>2}/MERRA2_{spinup}.tavg1_2d_flx_Nx.{year}{month:0>2}{day:0>2}.nc4',
		fn = os.path.join(merra2_dir, '{year}/{month:0>2}/MERRA2_{spinup}.tavg1_2d_flx_Nx.{year}{month:0>2}{day:0>2}.nc4'),
		variables = ['ustar','z0m','disph','rhoa','ulml','vlml','tstar','hlml','tlml','pblh','hflux','eflux']
	),
	'surface_flux_monthly': dict(
		file_granularity="monthly",
		tasks_func=tasks_monthly_merra2,
		meta_prepare_func=prepare_meta_merra2,
		prepare_func=prepare_month_surface_flux,
		template=os.path.join(merra2_dir, '{year}/MERRA2_*.tavgM_2d_flx_Nx.*.nc4'),
	    url = 'https://goldsmr4.gesdisc.eosdis.nasa.gov/data/MERRA2_MONTHLY/M2TMNXFLX.5.12.4/{year}/MERRA2_{spinup}.tavgM_2d_flx_Nx.{year}{month:0>2}.nc4',
		fn = os.path.join(merra2_dir, '{year}/MERRA2_{spinup}.tavgM_2d_flx_Nx.{year}{month:0>2}.nc4'),
		variables = ['ustar','z0m','disph','rhoa','ulml','vlml','tstar','hlml','tlml','pblh','hflux','eflux']
	),
	'surface_flux_dailymeans': dict(
		file_granularity="dailymeans",
		tasks_func=tasks_daily_merra2,
		meta_prepare_func=prepare_meta_merra2,
		prepare_func=prepare_dailymeans_surface_flux,
		template=os.path.join(merra2_dir, '{year}/{month:0>2}/MERRA2_*.statD_2d_slv_Nx.*.nc4'),
	    url = 'https://goldsmr4.gesdisc.eosdis.nasa.gov/data/MERRA2/M2SDNXSLV.5.12.4/{year}/{month:0>2}/MERRA2_{spinup}.statD_2d_slv_Nx.{year}{month:0>2}{day:0>2}.nc4',
		fn = os.path.join(merra2_dir, '{year}/{month:0>2}/MERRA2_{spinup}.statD_2d_slv_Nx.{year}{month:0>2}{day:0>2}.nc4'),
		variables = ["HOURNORAIN", "TPRECMAX","T2MMAX", "T2MMEAN", "T2MMIN"]
	)
}

# list of routines in weather_data_config to download wind data
wind_files = ['surface_flux']

# TODO: same for solar
solar_files = []

## Whatever is calling this needs to be directed to the correct weather config instead
#meta_data_config = dict(prepare_func=prepare_meta_merra2,
#						 template=os.path.join(merra2_dir, '{year}/{month:0>2}/MERRA2_*.tavg1_2d_flx_Nx.*.nc4'))

# Separate files for each day (coded in weather_data_config list)
#daily_files = True # needs to be specified somewhere else

# Latitude stored south to north (ie forward, = True) or north to south
lat_direction = True

# Spinup variable
spinup_var = True
def spinup_year(year):
	if (year>=1980 and year<1992):
		spinup = '100'
	elif (year>=1992 and year<2001):
		spinup = '200'
	elif (year>=2001 and year<2011):
		spinup = '300'
	elif (year>=2011):
		spinup = '400'
	return spinup
