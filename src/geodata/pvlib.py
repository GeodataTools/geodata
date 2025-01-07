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
import xarray as xr
from timezonefinder import TimezoneFinder
from pvlib import pvsystem
from pvlib.location import Location
from pvlib.modelchain import ModelChain
from pvlib.solarposition import get_solarposition
logger = logging.getLogger(__name__)

__all__ = ["_prepare_pvlib_ds"]

from .convert import (
    get_var,
    convert_temperature
)

class ModelChainConfig:
    """
    Defines pvlib ModelChain parameters as a class that
    can be passed to one or more instances of pvlib_model().
    Allows user to reuse a common set of ModelChain parameters across multiple
    PVSystems or even multiple cutouts.  

    Parameters
    ----------
    clearsky_model : string, default 'ineichen'
        Specifies the clear-sky model. Passed to location.get_clearsky. 
        Only used when DNI is not found in the weather inputs.
    transposition_model : string, default 'haydavies'
        Specifies the transposition model. Passed to system.get_irradiance.
    solar_position_method : string, default 'nrel_numpy'
        Specifies the method for calculating solar positions. Passed to location.get_solarposition.
    airmass_model : string, default 'kastenyoung1989'
        Specifies the airmass model. Passed to location.get_airmass.
    dc_model : string or function, optional
        Specifies the DC model. Valid strings are 'sapm', 'desoto', 'cec', 'pvsyst', 'pvwatts'. 
        If not specified, the model will be inferred from the parameters of system.arrays[i].module_parameters. 
        A user-defined function may also be provided, with the ModelChain instance passed as the first argument.
    ac_model : string or function, optional
        Specifies the AC model. Valid strings are 'sandia', 'adr', 'pvwatts'. 
        If not specified, the model will be inferred from the parameters of system.inverter_parameters. 
        A user-defined function may also be provided, with the ModelChain instance passed as the first argument.
    aoi_model : string or function, optional
        Specifies the angle of incidence (AOI) model. Valid strings are 'physical', 'ashrae', 'sapm', 'martin_ruiz', 
        'interp', 'no_loss'. If not specified, the model will be inferred from the parameters of 
        system.arrays[i].module_parameters. A user-defined function may also be provided, 
        with the ModelChain instance passed as the first argument.
    spectral_model : string or function, optional
        Specifies the spectral model. Valid strings are 'sapm', 'first_solar', 'no_loss'. 
        If not specified, the model will be inferred from the parameters of system.arrays[i].module_parameters. 
        A user-defined function may also be provided, with the ModelChain instance passed as the first argument.
    temperature_model : string or function, optional
        Specifies the temperature model. Valid strings are 'sapm', 'pvsyst', 'faiman', 'fuentes', 'noct_sam'. 
        A user-defined function may also be provided, with the ModelChain instance passed as the first argument.
    dc_ohmic_model : string or function, default 'no_loss'
        Specifies the DC ohmic loss model. Valid strings are 'dc_ohms_from_percent', 'no_loss'. 
        A user-defined function may also be provided, with the ModelChain instance passed as the first argument.
    losses_model : string or function, default 'no_loss'
        Specifies the losses model. Valid strings are 'pvwatts', 'no_loss'. 
        A user-defined function may also be provided, with the ModelChain instance passed as the first argument.
    name : string, optional
        Specifies the name of the ModelChain instance.
    
    See also:
    - pvlib.modelchain.ModelChain(): 
        https://pvlib-python.readthedocs.io/en/stable/reference/generated/pvlib.modelchain.ModelChain.html

    """
    def __init__(
        self,
        clearsky_model='ineichen',
        transposition_model='haydavies',
        solar_position_method='nrel_numpy',
        airmass_model='kastenyoung1989',
        dc_model=None,
        ac_model=None,
        aoi_model=None,
        spectral_model=None,
        temperature_model=None,
        dc_ohmic_model='no_loss',
        losses_model='no_loss',
        name=None
    ):
        self.clearsky_model = clearsky_model
        self.transposition_model = transposition_model
        self.solar_position_method = solar_position_method
        self.airmass_model = airmass_model
        self.dc_model = dc_model
        self.ac_model = ac_model
        self.aoi_model = aoi_model
        self.spectral_model = spectral_model
        self.temperature_model = temperature_model
        self.dc_ohmic_model = dc_ohmic_model
        self.losses_model = losses_model
        self.name = name

    def model_chain_to_kwargs(self):
        return self.__dict__

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

def pv_system(*args, **kwargs):
    """
    Wrapper for pvlib.pvsystem.PVSystem().
    The PVSystem class defines a standard set of PV system attributes
    and modeling functions. This class describes the collection and 
    interactions of PV system components rather than an installed system
    on the ground. It is typically used in combination with Location 
    and ModelChain objects.

    The class supports basic system topologies consisting of:
        - N total modules arranged in series (modules_per_string=N, strings_per_inverter=1).
        - M total modules arranged in parallel (modules_per_string=1, strings_per_inverter=M).
        - NxM total modules arranged in M strings of N modules each 
          (modules_per_string=N, strings_per_inverter=M).

    Parameters
    ----------
    arrays : array (optional)
            An Array or list of arrays that are part of the system. 
            See pvlib documentation for full description.
    surface_tilt : float
            Surface tilt angles in decimal degrees. The tilt angle is 
            defined as degrees from horizontal (e.g. surface facing up = 0, 
            surface facing horizon = 90).
    surface_azimuth : float
            Azimuth angle of the module surface. North=0, East=90, South=180, West=270.
    albedo : float
            Ground surface albedo. If not supplied, then surface_type is used to look up 
            a value in pvlib.albedo.SURFACE_ALBEDOS. If surface_type is also not supplied 
            then a ground surface albedo of 0.25 is used.
    surface_type : string
            The ground surface type. See pvlib.albedo.SURFACE_ALBEDOS for valid values.
    module : string
            The model name of the modules. May be used to look up the module_parameters dictionary via some other method.
    module_type : string 
            Describes the module’s construction. Valid strings are ‘glass_polymer’ and ‘glass_glass’. 
            Used for cell and module temperature calculations.
    module_parameters : dict
            Module parameters as defined by the SAPM, CEC, or other.
    temperature_model_parameters : dict
            Temperature model parameters as required by one of the models in pvlib.temperature (excluding poa_global, temp_air and wind_speed).
    modules_per_string : int, float
            See system topology discussion above.
    strings_per_inverter : int, float 
            See system topology discussion above.
    inverter : string 
            The model name of the inverters. May be used to look up the inverter_parameters dictionary via some other method.
    inverter_parameters : dict
            Inverter parameters as defined by the SAPM, CEC, or other.
    racking_model : string 
            Valid strings are ‘open_rack’, ‘close_mount’, and ‘insulated_back’. 
            Used to identify a parameter set for the SAPM cell temperature model.
    losses_parameters : dict 
            Losses parameters as defined by PVWatts or other.    
    name : string (optional)
    
    See also: https://pvlib-python.readthedocs.io/en/stable/reference/generated/pvlib.pvsystem.PVSystem.html
    """
    return pvsystem.PVSystem(*args, **kwargs)

def pvlib_model(
        cutout, 
        system,
        model_chain_config,
        vars = [
            'influx_diffuse', 
            'influx_direct', 
            'temperature', 
            'wnd100m',
            'total_column_water_vapour'
        ]
    ):

    """
    Wrapper function that applies a pvlib model as defined by pvlib.modelchain.ModelChain across all unique coordinates
    represented by a geodata cutout.  

    Parameters
    ----------

    cutout : cutout
            Geodata cutout based on ERA5 dataset.
    system : system

    """

    weather_data = _prepare_pvlib_ds(cutout, *vars).to_dataframe()
    unique_coords = weather_data.index.droplevel('time').drop_duplicates()
    coord_subsets = []
    for y, x in unique_coords:
        subset = weather_data.loc[(slice(None), y, x), :].reset_index(['x', 'y'])
        tz_str = TimezoneFinder().timezone_at(lat=y, lng=x)
        location = Location(latitude=y, longitude=x, tz = tz_str)
        
        mc = ModelChain(
            system, 
            location, 
            **model_chain_config.model_chain_to_kwargs()
        )
        mc.run_model(subset)
        
        subset['ac'] = mc.results.ac
        subset.loc[subset['ac'] < 0, 'ac'] = 0
        subset = subset.merge(mc.results.dc, how = "left", on = "time").set_index(['x', 'y'], append=True)

        coord_subsets.append(subset)

    weather_data_final = pd.concat(coord_subsets)

    return xr.Dataset.from_dataframe(weather_data_final)



## solar PV - pvlib
def _calculate_pvlib_solarposition(ds):
    nt, ny, nx = ds.sizes['time'], ds.sizes['y'], ds.sizes['x']
    time_expanded = np.broadcast_to(ds.time.values[:, None, None], (nt, ny, nx)).ravel()
    yy, xx = np.meshgrid(ds.y, ds.x, indexing="ij")
    x_expanded = np.tile(xx.ravel(), nt)
    y_expanded = np.tile(yy.ravel(), nt)
    solarposition = get_solarposition(time_expanded, y_expanded, x_expanded)
    solarposition.index = pd.MultiIndex.from_arrays([time_expanded, y_expanded, x_expanded], names=['time', 'y', 'x'])
    return solarposition

def _calculate_ghi(ds, zenith):
    dhi = ds.influx_diffuse.values.ravel()
    dni = ds.influx_direct.values.ravel()
    ghi = np.clip(
        dhi + dni * np.cos(zenith),
        0,
        np.Inf
    )

    reshaped_ghi = ghi.values.reshape(
        ds.sizes['time'], 
        ds.sizes['y'], 
        ds.sizes['x']
    )
    
    ghi = xr.DataArray(
        reshaped_ghi,
        dims=("time", "y", "x"),
        coords={
            "time": ds['time'].values, 
            "y": ds['y'].values, 
            "x": ds['x'].values
        },
        name="ghi"
    )
    return ghi


def _prepare_pvlib_ds(cutout, *varnames):
    ds = xr.Dataset({
        name: get_var(cutout, name)
        for name in varnames
    })
    
    sp = _calculate_pvlib_solarposition(ds)
    ghi = _calculate_ghi(ds, sp['zenith'])

    ds = (
        ds
        .assign(
            ghi=ghi,
            temperature=convert_temperature(ds),
        )
        .rename({
            'influx_diffuse': 'dhi',
            'influx_direct': 'dni',
            'temperature': 'temp_air',
            'wnd100m': 'wind_speed',
            'total_column_water_vapour': 'precipitable_water'
        })
    )

    return ds

