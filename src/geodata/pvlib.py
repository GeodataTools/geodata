# Copyright 2020 Michael Davidson (UCSD).

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


"""
GEODATA

Geospatial Data Collection and "Pre-Analysis" Tools
"""

import logging
import pandas as pd
import numpy as np
from pvlib.solarposition import get_solarposition
from pvlib import pvsystem
logger = logging.getLogger(__name__)

__all__ = ["_prepare_pvlib_df"]

from .convert import (
    get_var,
    convert_temperature
)

def retrieve_sam(samfile, path=None):
    """
    Wrapper for pvlib.pvsystem.retrieve_sam(). Retrieves latest module 
    and inverter info from a file bundled with pvlib, a path or a 
    URL (like SAM’s website), and returns it as a Pandas DataFrame.

    Supported databases:
    - CEC module database
    - Sandia Module database
    - CEC Inverter database
    - Anton Driesse Inverter database

    Parameters
    ----------
    name : string
            Use one of the following strings to retrieve a database bundled with pvlib:
                - ’CECMod’ - returns the CEC module database
                - ’CECInverter’ - returns the CEC Inverter database
                - ’SandiaInverter’ - returns the CEC Inverter database 
                    (CEC is only current inverter db available; tag kept for backwards compatibility)
                - ’SandiaMod’ - returns the Sandia Module database
                - ’ADRInverter’ - returns the ADR Inverter database
    
    Optional Parameters
    ----------
    path : string
            Path to a CSV file or a URL.

    Returns: DataFrame

    See also:
        - pvlib.pvsystem.retrieve_sam(): 
            https://pvlib-python.readthedocs.io/en/stable/reference/generated/pvlib.pvsystem.retrieve_sam.html

    """
    return pvsystem.retrieve_sam(name=samfile, path=path)

## solar PV - pvlib
def _calculate_pvlib_solarposition(time, y, x):
    sp = get_solarposition(time, y, x)
    sp.index = pd.MultiIndex.from_arrays([time, y, x], names=['time', 'y', 'x'])
    return sp

def _calculate_ghi(dhi, dni, zenith):
    return np.clip(
        dhi + dni * np.cos(zenith),
        0,
        np.Inf
    )

def _prepare_pvlib_df(cutout, *args):
    dfs = []
    for varname in args:
        
        var_array =  get_var(cutout, varname)
        df = var_array.to_dataframe(name=var_array.name)
        dfs.append(df)
    
    result_df = pd.concat(dfs, axis=1)
    result_df.reset_index(inplace=True)
    result_df.set_index(['time', 'y', 'x'], inplace=True)
    result_df = (
        result_df
        .loc[:, ~result_df.columns.duplicated()]
        .assign(temperature=convert_temperature(result_df))  # Convert temperature to C
        .rename(columns={
            'influx_diffuse': 'dhi',
            'influx_direct': 'dni',
            'temperature': 'temp_air',
            'wnd100m': 'wind_speed'
        })
    )    
    # solar position
    sp = _calculate_pvlib_solarposition(
        result_df.index.get_level_values('time'),
        result_df.index.get_level_values('y'),
        result_df.index.get_level_values('x')
    )

    # ghi
    result_df['ghi'] = _calculate_ghi(
        result_df['dhi'], 
        result_df['dni'], 
        sp['zenith']
    )

    # precipitation
    result_df['precipitable_water'] = 0.5 # guesstimate from example script, should find a better way to call this
    return result_df