# -*- coding: utf-8 -*-
## wind.py
#	functions to extrapolate wind speeds from certain heights to other heights

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
Renewable Energy Atlas Lite (Atlite)

Light-weight version of Aarhus RE Atlas for converting weather data to power systems data
"""
import xarray as xr
import numpy as np

import logging
logger = logging.getLogger(__name__)


"""
Extrapolation functions
	Note: these are called internally
"""

def log_ratio(ds, to_height, from_height, from_name):
	""" Logarithmic ratio law
	# Equation (2) in Andresen, G. et al (2015):
	# 	'Validation of Danish wind time series from a new global renewable
	# 	energy atlas for energy system analysis'.
	"""

	wnd_spd = ds[from_name] * ( np.log(to_height /ds['roughness'])
							  / np.log(ds[from_height]/ds['roughness']))

	wnd_spd.attrs.update({"long_name":
							f"extrapolated {to_height} m wind speed using log ratio",
						  "units" : "m s**-1"})
	return wnd_spd

def log_law(ds, to_height, from_height, from_name):
	""" Logarithmic (integration) law
		[3] S. Emeis, Wind Energy Meteorology (Springer, Berlin, 2013).
	"""
	vonk = 0.4
	wnd_spd = ds[from_name] + 	( ds['ustar'] / vonk
								* np.log((to_height - ds['disph']) / ds[from_height]) )

	wnd_spd.attrs.update({"long_name":
							f"extrapolated {to_height} m wind speed using log (integration) law",
						  "units" : "m s**-1"})
	return wnd_spd


"""
Main call (from convert.convert_wind)
"""

def extrapolate_wind_speed(ds, to_height, extrap_fn = log_ratio, from_height=None, var_height=None):
	"""Extrapolate the wind speed from a given height above ground to another.

	If ds already contains a key refering to wind speeds at the desired to_height,
	no conversion is done and the wind speeds are directly returned.

	Otherwise, extrapolates according to (1) extrap_fn and (2) heights



	Parameters
	----------
	ds : xarray.Dataset
		Dataset containing the wind speed time-series
	to_height : int|float
		Height (m) to which the wind speeds are extrapolated
	extrap_fn : function for wind speed extrapolation
		log_ratio : wind speed follows the logarithmic ratio law as desribed in [1]
		log_law : wind speed follows logarithmic (integration) law described in [3]
		power_law : wind speed follows power law (with fixed alpha), e.g., in [4]
	from_height : int
		(Optional)
		Height (m) from which the wind speeds are interpolated to 'to_height'.
		If not provided, the closest height to 'to_height' is selected.
	var_height : str
		(Optional)
		suffix of variables in ds corresponding to variable height
		e.g., `lml` => height contained in `hlml`, wind speed contained in `wndlml`

	Returns
	-------
	da : xarray.DataArray
		DataArray containing the extrapolated wind speeds. Name of the DataArray
		is 'wnd{to_height:d}'.

	References
	----------
	[1] Equation (2) in Andresen, G. et al (2015):
		'Validation of Danish wind time series from a new global renewable
		energy atlas for energy system analysis'.
	[2] https://en.wikipedia.org/w/index.php?title=Roughness_length&oldid=862127433,
		Retrieved 2019-02-15.
	[3] S. Emeis, Wind Energy Meteorology (Springer, Berlin, 2013).
	[4] Archer, C.L., Jacobson, M.Z., 2005. Evaluation of global wind power.
		Journal of Geophysical Research 110, D12110.
	"""


	to_name   = "wnd{h:0d}m".format(h=int(to_height))
	if to_name in ds:
		# already found wind speed at given height in dataset
		return ds[to_name]

	# Sanitize roughness for logarithm: 0.0002 corresponds to open water [2]
	ds['roughness'].values[ds['roughness'].values <= 0.0] = 0.0002

	if not from_height is None:
		# passed a from_height
		if not var_height is None:
			raise AssertionError("Cannot pass both from_height and var_height to extrapolate_wind_speed")
		from_name = "wnd{h:0d}m".format(h=int(from_height))
		ds['from_height'] = from_height

		wnd_spd = extrap_fn(ds, to_height, 'from_height', from_name)
		wnd_spd.attrs['long_name'] = wnd_spd.attrs['long_name'] + ', ' + f'from fixed height = {from_height}'

	elif not var_height is None:
		# passed a variable height (eg lml)
		# set variable names
		from_height = f"h{var_height}"
		from_name = f"wnd{var_height}"

		wnd_spd = extrap_fn(ds, to_height, from_height, from_name)
		wnd_spd.attrs['long_name'] = wnd_spd.attrs['long_name'] + ', ' + f'from variable height = {var_height}'

	else:
		# based on nearest height
		heights = np.asarray([int(s[3:-1]) for s in ds if s.startswith("wnd")])
		if len(heights) == 0:
			raise AssertionError("Wind speed is not in dataset")

		from_height = heights[np.argmin(np.abs(heights-to_height))]
		from_name = "wnd{h:0d}m".format(h=int(from_height))
		ds['from_height'] = from_height

		wnd_spd = extrap_fn(ds, to_height, 'from_height', from_name)
		wnd_spd.attrs['long_name'] = wnd_spd.attrs['long_name'] + ', ' + f'from nearest height = {from_height}'

	return wnd_spd.rename(to_name)
