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
import logging
import numpy as np
import xarray as xr

from ..config import era5_dir

logger = logging.getLogger(__name__)
datadir = era5_dir

try:
	import cdsapi
	has_cdsapi = True
except ImportError:
	has_cdsapi = False

# Model and Projection Settings
projection = 'latlong'

def _rename_and_clean_coords(ds, add_lon_lat=True):
	"""Rename 'lon'/'longitude' and 'lat'/'latitude' columns to 'x' and 'y'

	Optionally (add_lon_lat, default:True) preserves latitude and longitude columns as 'lat' and 'lon'.
	"""
	# Rename latitude / lat -> y, longitude / lon -> x
	if 'latitude' in list(ds.coords):
		ds = ds.rename({'latitude': 'y'})
	if 'longitude' in list(ds.coords):
		ds = ds.rename({'longitude': 'x'})
	if 'lat' in list(ds.coords):
		ds = ds.rename({'lat': 'y'})
	if 'lon' in list(ds.coords):
		ds = ds.rename({'lon': 'x'})

	if add_lon_lat:
		ds = ds.assign_coords(lon=ds.coords['x'], lat=ds.coords['y'])
	return ds

def api_hourly_era5(
	toDownload,
	bounds,
	download_vars,
	product,
	product_type
	):
	if not has_cdsapi:
		raise RuntimeError(
					"Need installed cdsapi python package available from "
					"https://cds.climate.copernicus.eu/api-how-to"
				)

	if len(toDownload) == 0:
		logger.info("All ERA5 files for this dataset have been downloaded.")
	else:
		logger.info("Preparing to download %s files.", str(len(toDownload)))

		for f in toDownload:
			print(f)
			os.makedirs(os.path.dirname(f[1]), exist_ok=True)

			## for each file in self.todownload - need to then reextract year month in order to make query
			query_year = str(f[2])
			query_month = str(f[3]) if len(str(f[3])) == 2 else '0' + str(f[3])

			#2. Full data file
			full_request = {
				'product_type': product_type,
				'format':'netcdf',
				'year':query_year,
				'month':query_month,
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
				],
				'variable': download_vars
			}

			if bounds is not None:
				full_request['area'] = bounds

			full_result = cdsapi.Client().retrieve(
				product,
				full_request
			)

			logger.info("Downloading metadata request for %s variables to %s", len(full_request['variable']), f)
			full_result.download(f[1])
			logger.info("Successfully downloaded to %s", f[1])

def api_monthly_era5(
	toDownload,
	bounds,
	download_vars,
	product,
	product_type
	):
	if not has_cdsapi:
		raise RuntimeError(
					"Need installed cdsapi python package available from "
					"https://cds.climate.copernicus.eu/api-how-to"
				)

	if len(toDownload) == 0:
		logger.info("All ERA5 files for this dataset have been downloaded.")
	else:
		logger.info("Preparing to download %s files.", str(len(toDownload)))

		for f in toDownload:
			print(f)
			os.makedirs(os.path.dirname(f[1]), exist_ok=True)

			## for each file in self.todownload - need to then reextract year month in order to make query
			query_year = str(f[2])
			query_month = str(f[3]) if len(str(f[3])) == 2 else '0' + str(f[3])

			#2. Full data file
			full_request = {
				'product_type':product_type,
				'format':'netcdf',
				'year':query_year,
				'month':query_month,
				'time':'00:00',
				'variable': download_vars
			}

			if bounds is not None:
				full_request['area'] = bounds

			full_result = cdsapi.Client().retrieve(
				product,
				full_request
			)

			logger.info("Downloading metadata request for %s variables to %s", len(full_request['variable']), f)
			full_result.download(f[1])
			logger.info("Successfully downloaded to %s", f[1])



def convert_and_subset_lons_lats_era5(ds, xs, ys):
	# Rename geographic dimensions to x,y
	# Subset x,y according to xs, ys (subset_x_y_era5)

	# Rename lat and lon
	ds = _rename_and_clean_coords(ds)

	# Longitudes should go from -180. to +180.
	if len(ds.coords['x'].sel(x=slice(xs.start + 360., xs.stop + 360.))):
		ds = xr.concat([ds.sel(x=slice(xs.start + 360., xs.stop + 360.)),
						ds.sel(x=xs)],
					   dim="x")
		ds = ds.assign_coords(x=np.where(ds.coords['x'].values <= 180,
											 ds.coords['x'].values,
											 ds.coords['x'].values - 360.))
	# Subset x and y
	ds = subset_x_y_era5(ds, xs, ys)

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

def prepare_meta_era5(xs, ys, year, month, template):
	# Reference of the quantities
	# https://confluence.ecmwf.int/display/CKB/ERA5+data+documentation
	# Geopotential is aka Orography in the CDS:
	# https://confluence.ecmwf.int/pages/viewpage.action?pageId=78296105

	fns = glob.iglob(template.format(year=year, month=month))
	try:
		with xr.open_mfdataset(fns, combine='by_coords') as ds0:
			ds = ds0.coords.to_dataset()
			ds = convert_and_subset_lons_lats_era5(ds, xs, ys)
			meta = ds.load()
	except Exception as e:
		logger.exception("Error when preparing for cutout: %s",
						 e.args[0])
		raise e
	return meta

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
		logger.info('Opening %s', fn)
		ds = _rename_and_clean_coords(ds)
		# ds = _add_height(ds)
		ds = subset_x_y_era5(ds, xs, ys)

		ds = ds.rename({'fdir': 'influx_direct', 'tisr': 'influx_toa'})
		with np.errstate(divide='ignore', invalid='ignore'):
			ds['albedo'] = (((ds['ssrd'] - ds['ssr'])/ds['ssrd']).fillna(0.)
							.assign_attrs(units='(0 - 1)', long_name='Albedo'))
		ds['influx_diffuse'] = ((ds['ssrd'] - ds['influx_direct'])
								.assign_attrs(units='J m**-2',
											long_name='Surface diffuse solar radiation downwards'))
		ds = ds.drop(['ssrd', 'ssr'])

		# Convert from energy to power J m**-2 -> W m**-2 and clip negative fluxes
		for a in ('influx_direct', 'influx_diffuse', 'influx_toa'):
			ds[a] = ds[a].clip(min=0.) / (60.*60.)
			ds[a].attrs['units'] = 'W m**-2'

		ds['wnd100m'] = (np.sqrt(ds['u100']**2 + ds['v100']**2)
						.assign_attrs(units=ds['u100'].attrs['units'],
									long_name="100 metre wind speed"))
		ds = ds.drop(['u100', 'v100'])

		ds = ds.rename({'ro': 'runoff',
						't2m': 'temperature',
						'sp': 'pressure',
						'stl4': 'soil temperature',
						'fsr': 'roughness'
						})

		ds['runoff'] = ds['runoff'].clip(min=0.)

		yield (year, month), ds


def tasks_monthly_era5(xs, ys, yearmonths, prepare_func, **meta_attrs):
	if not isinstance(xs, slice):
		xs = slice(*xs.values[[0, -1]])
	if not isinstance(ys, slice):
		ys = slice(*ys.values[[0, -1]])
	fn = meta_attrs['fn']

	logger.info(yearmonths)
	logger.info(list(yearmonths))

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
	'wind_solar_hourly': dict(
		api_func=api_hourly_era5,
		file_granularity="monthly",
		tasks_func=tasks_monthly_era5,
		meta_prepare_func=prepare_meta_era5,
		prepare_func=prepare_month_era5,
		template=os.path.join(era5_dir, '{year}/{month:0>2}/wind_solar_hourly.nc'),
		fn = os.path.join(era5_dir, '{year}/{month:0>2}/wind_solar_hourly.nc'),
		product='reanalysis-era5-single-levels',
		product_type='reanalysis',
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
					   'total_sky_direct_solar_radiation_at_surface',
					   'forecast_surface_roughness',
					   'orography'
				   ]
		),
	'wind_solar_monthly': dict(
		api_func=api_monthly_era5,
		file_granularity="monthly",
		tasks_func=tasks_monthly_era5,
		meta_prepare_func=prepare_meta_era5,
		prepare_func=prepare_month_era5,
		template=os.path.join(era5_dir, '{year}/{month:0>2}/wind_solar_monthly.nc'),
		fn = os.path.join(era5_dir, '{year}/{month:0>2}/wind_solar_monthly.nc'),
		product='reanalysis-era5-single-levels-monthly-means',
		product_type='monthly_averaged_reanalysis',
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
					   'total_sky_direct_solar_radiation_at_surface',
					   'forecast_surface_roughness',
					   'orography'
				   ]
		)
}

# No separate files for each day (would be coded in weather_data_config list, see merra2.py)
daily_files = False

# Latitude direction stored
# 	South to north = True
#	North to south = False
lat_direction = False

# Spinup variable (necessary for MERRA)
spinup_var = False
