# Copyright 2016-2017 Gorm Andresen (Aarhus University), Jonas Hoersch (FIAS), Tom Brown (FIAS)

# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License as
# published by the Free Software Foundation; either version 3 of the
# License, or (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


"""
This module contains various functions used to perform conversion in Geodata.
"""

import logging
from operator import itemgetter
from typing import TypeVar

import numpy as np
import xarray as xr

from . import wind as windm
from .pv.irradiation import TiltedIrradiation
from .pv.orientation import SurfaceOrientation, get_orientation  # noqa: F401
from .pv.solar_panel_model import SolarPanelModel
from .pv.solar_position import SolarPosition
from .resource import get_solarpanelconfig, get_windturbineconfig, windturbine_smooth  # noqa: F401

logger = logging.getLogger(__name__)
Cutout = TypeVar("Cutout")


def convert_solar_thermal(
    ds, orientation, trigon_model, clearsky_model, c0, c1, t_store
):
    # convert storage temperature to Kelvin in line with reanalysis data
    t_store += 273.15

    # Downward shortwave radiation flux is in W/m^2
    # http://rda.ucar.edu/datasets/ds094.0/#metadata/detailed.html?_do=y
    solar_position = SolarPosition(ds)
    surface_orientation = SurfaceOrientation(ds, solar_position, orientation)
    irradiation = TiltedIrradiation(
        ds, solar_position, surface_orientation, trigon_model, clearsky_model
    )

    # overall efficiency; can be negative, so need to remove negative values below
    eta = c0 - c1 * ((t_store - ds["temperature"]) / irradiation)

    output = irradiation * eta

    return (output).where(output > 0.0).fillna(0.0)


def convert_pv(ds, panel, orientation, trigon_model="simple", clearsky_model="simple"):
    solar_position = SolarPosition(ds)
    surface_orientation = SurfaceOrientation(ds, solar_position, orientation)
    irradiation = TiltedIrradiation(
        ds,
        solar_position,
        surface_orientation,
        trigon_model=trigon_model,
        clearsky_model=clearsky_model,
    )
    solar_panel = SolarPanelModel(ds, irradiation, panel)
    return solar_panel


## wind


def convert_wind(ds, turbine, **params):
    """
    Convert wind speeds for turbine to wind energy generation.
    Selects hub height according to turbine model

    - load turbine parameters
    - extrapolate wind speeds 			(wind.extrapolate_wind_speed)
            extrapolate_wind_speed(ds, to_height, extrap_fn = log_ratio, from_height=None, var_height=None)

    Optional Parameters
    ------

    extrap_fn : function for extrapolation
    from_height (int) : fixed height from which to extrapolate
    var_height (str) : suffix for variables containing wind speed and variable height

    """

    V, POW, hub_height, P = itemgetter("V", "POW", "hub_height", "P")(turbine)

    wnd_hub = windm.extrapolate_wind_speed(ds, to_height=hub_height, **params)

    return xr.DataArray(np.interp(wnd_hub, V, POW / P), coords=wnd_hub.coords)


def convert_windspd(ds, hub_height, **params):
    """
    Extract wind speeds at given height

    - extrapolate wind speeds 			(wind.extrapolate_wind_speed)
            extrapolate_wind_speed(ds, to_height, extrap_fn = log_ratio, from_height=None, var_height=None)

    Parameters
    ----------
    hub_height : num
            extrapolation height

    Optional Parameters
    ------

    extrap_fn : function for extrapolation
    from_height (int) : fixed height from which to extrapolate
    var_height (str) : suffix for variables containing wind speed and variable height

    """
    wnd_hub = windm.extrapolate_wind_speed(ds, to_height=hub_height, **params)

    return xr.DataArray(wnd_hub, coords=wnd_hub.coords)


def convert_windwpd(ds, hub_height, **params):
    """
    Extract wind power density at given height, according to:
            WPD = 0.5 * Density * Windspd^3

    - extrapolate wind speeds 			(wind.extrapolate_wind_speed)
            extrapolate_wind_speed(ds, to_height, extrap_fn = log_ratio, from_height=None, var_height=None)

    Parameters
    ----------
    hub_height : num
            extrapolation height

    Optional Parameters
    ------

    extrap_fn : function for extrapolation
    from_height (int) : fixed height from which to extrapolate
    var_height (str) : suffix for variables containing wind speed and variable height

    """
    wnd_hub = windm.extrapolate_wind_speed(ds, to_height=hub_height, **params)

    return xr.DataArray(0.5 * ds["rhoa"] * wnd_hub**3, coords=wnd_hub.coords)


def convert_pm25(ds):
    """
    Generate PM2.5 time series according to [1]:

            PM2.5 = [Dust2.5] + [SS2.5] + [BC] + 1.4*[OC] + 1.375*[SO4]

    Parameters
    ----------
    **params : None needed currently.

    References
    -------
    [1] Buchard, V., da Silva, A. M., Randles, C. A., Colarco, P., Ferrare, R., Hair, J., … Winker, D. (2016).
        Evaluation of the surface PM2.5 in Version 1 of the NASA MERRA Aerosol Reanalysis
        over the United States. Atmospheric Environment, 125, 100-111.
    https://doi.org/10.1016/j.atmosenv.2015.11.004
    """

    ds["pm25"] = (
        ds["dusmass25"]
        + ds["sssmass25"]
        + ds["bcsmass"]
        + 1.4 * ds["ocsmass"]
        + 1.375 * ds["so4smass"]
    )

    return 1e9 * ds["pm25"]  # kg / m3 to ug / m3


# Manipulate arbitrary variables


def _get_var(ds, var):
    """
    (Internal) Extract a specific variable from cutout
    See: get_var
    """
    return xr.DataArray(ds[var], coords=ds.coords)


def get_var(cutout, var, **params):
    """
    Extract a specific variable from cutout

    Parameters
    ----------
    var : str
            Name of variable to extract from dataset

    Returns: dataarray
    """
    logger.info("Getting variable: %s", str(var))
    return cutout.convert_cutout(convert_func=_get_var, var=var, **params)


def _compute_var(ds, fn):
    """
    (Internal) Compute a specific function from cutout
    See: compute_var
    """
    return xr.DataArray(fn(ds), coords=ds.coords)


def compute_var(cutout, fn, **params):
    """
    Compute a specific function from cutout

    Parameters
    ----------
    var : str
            Name of variable to extract from dataset

    Returns: dataarray
    """
    logger.info("Computing variable: %s", str(fn))
    return cutout.convert_cutout(convert_func=_compute_var, fn=fn, **params)
