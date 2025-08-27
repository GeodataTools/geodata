# Copyright 2016-2017 Gorm Andresen (Aarhus University), Jonas Hoersch (FIAS), Tom Brown (FIAS)
# Copyright 2025 Xiqiang Liu, Michael Davidson (UCSD)

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
This module contains various functions used to perform conversion in Geodata.
"""

import datetime as dt
import logging
from operator import itemgetter
from typing import TYPE_CHECKING, Literal, Callable

import numpy as np
import xarray as xr

from . import wind as windm
from .pv.irradiation import TiltedIrradiation
from .pv.orientation import SurfaceOrientation, get_orientation  # noqa: F401
from .pv.solar_panel_model import SolarPanelModel
from .pv.solar_position import SolarPosition
from .resource import (
    get_solarpanelconfig,
    get_windturbineconfig,
    windturbine_smooth,
)

if TYPE_CHECKING:
    from .cutout import Cutout
else:
    Cutout = object

logger = logging.getLogger(__name__)


# Heat Demand
def convert_heat_demand(
    ds: xr.Dataset,
    threshold: float,
    a: float,
    constant: float,
    hour_shift: float,
):
    # Temperature is in Kelvin; take daily average
    T = ds["temperature"]
    T.coords["time"] += np.timedelta64(dt.timedelta(hours=hour_shift))

    T = ds["temperature"].resample(time="1D").mean(dim="time")
    threshold += 273.15
    heat_demand_value = a * (threshold - T)

    heat_demand_value.values[heat_demand_value.values < 0.0] = 0.0

    return constant + heat_demand_value


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

    return output.where(output > 0.0).fillna(0.0)


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


def heat_demand(
    cutout: Cutout,
    threshold: float = 15.0,
    a: float = 1.0,
    constant: float = 0.0,
    hour_shift: float = 0.0,
    **params,
):
    """Convert outside temperature into daily heat demand using the
    degree-day approximation.

    Since "daily average temperature" means different things in
    different time zones and since xarray coordinates do not handle
    time zones gracefully like pd.DateTimeIndex, you can provide an
    hour_shift to redefine when the day starts.

    E.g. for Moscow in winter, hour_shift = 4, for New York in winter,
    hour_shift = -5

    This time shift applies across the entire spatial scope of ds for
    all times. More fine-grained control will be built in a some
    point, i.e. space- and time-dependent time zones.

    WARNING: Because the original data is provided every month, at the
    month boundaries there is untidiness if you use a time shift. The
    resulting xarray will have duplicates in the index for the parts
    of the day in each month at the boundary. You will have to
    re-average these based on the number of hours in each month for
    the duplicated day.

    Args:
        threshold (float): Outside temperature in degrees Celsius above which there is no heat demand.
        a (float): Linear factor relating heat demand to outside temperature.
        constant (float): Constant part of heat demand that does not depend on outside
            temperature (e.g. due to water heating).
        hour_shift (float): Time shift relative to UTC for taking daily average

    Returns:
        xr.DataArray: Heat demand

    Note:
        You can also specify all of the general conversion arguments
        documented in the `convert_cutout` function.
    """

    return cutout._convert_cutout(
        convert_func=convert_heat_demand,
        threshold=threshold,
        a=a,
        constant=constant,
        hour_shift=hour_shift,
        **params,
    )


def temperature(cutout: Cutout, **convert_params):
    """Convert temperature in Cutout to outside temperature.

    Args:
        convert_params: Keyword arguments passed to `convert_cutout` function

    Returns:
        xr.DataArray: Data of the Cutout with temperature converted to outside temperatures.
    """
    return cutout._convert_cutout(
        convert_func=lambda ds: ds["temperature"] - 273.15, **convert_params
    )


def soil_temperature(cutout: Cutout, **convert_params):
    """Return soil temperature (useful for e.g. heat pump T-dependent
    coefficient of performance).

    Args:
        convert_params: Keyword arguments passed to `convert_cutout` function

    Returns:
        xr.DataArray: Data of the Cutout with temperature converted to soil temperatures.
    """
    return cutout._convert_cutout(
        convert_func=lambda ds: (ds["soil temperature"] - 273.15).fillna(0.0),
        **convert_params,
    )


def solar_thermal(
    cutout: Cutout,
    orientation: dict | str | Callable | None = None,
    trigon_model: str = "simple",
    clearsky_model: Literal["simple", "enhanced"] = "simple",
    c0: float = 0.8,
    c1: float = 3.0,
    t_store: float = 80.0,
    **params,
):
    """Convert downward short-wave radiation flux and outside temperature
    into time series for solar thermal collectors.

    Mathematical model and defaults for c0, c1 based on model in [1].

    Args:
        orientation (Union[dict, str, callable]): Panel orientation with slope and azimuth
            (units of degrees), or 'latitude_optimal'.
        trigon_model (str): Type of trigonometry model
        clearsky_model (str): Type of clearsky model for diffuse irradiation. Either
            `simple` or `enhanced`.
        c0 (float): Parameter for model in [1] This defaults to 0.8.
        c1 (float): Parameter for model in [1] This defaults to 3.0.
        t_store (float): Store temperature in degree Celsius

    Note:
        You can also specify all of the general conversion arguments
        documented in the `convert_cutout` function.

    References:
        [1] Henning and Palzer, Renewable and Sustainable Energy Reviews 30
        (2014) 1003-1018
    """

    if orientation is None:
        orientation = {"slope": 45.0, "azimuth": 180.0}

    if not callable(orientation):
        orientation = get_orientation(orientation)

    return cutout._convert_cutout(
        convert_func=convert_solar_thermal,
        orientation=orientation,
        trigon_model=trigon_model,
        clearsky_model=clearsky_model,
        c0=c0,
        c1=c1,
        t_store=t_store,
        **params,
    )


def wind(
    cutout: Cutout,
    turbine: str | dict,
    method: Literal["simple", "interpolation", "extrapolation"],
    smooth: bool | dict = False,
    **params,
):
    """Convert wind speed time-series into wind generation time-series.

    Args:
        turbine (Union[str, dict]): Name of a turbine or a dictionary with the parameters
            for the wind turbine in [2].
        smooth (Union[bool, dict]): If True, the wind speed time-series will be smoothed
            before conversion. If False, no smoothing will be applied. If a dictionary is
            passed, the smoothing parameters will be used.
        **params: Keyword arguments passed to `convert_cutout` function
    """

    if isinstance(turbine, str):
        turbine = get_windturbineconfig(turbine)

    if smooth:
        turbine = windturbine_smooth(turbine, params=smooth)

    match method:
        case "simple":
            return cutout._convert_cutout(
                convert_func=convert_wind, turbine=turbine, **params
            )

        case _:
            raise ValueError(f"Method {method} not supported.")


def windspd(cutout: Cutout, **params):
    """
    Generate wind speed time-series

    convert.convert_cutout → convert.convert_windspd

    Parameters
    ----------
    **params
        Must have 1 of:
            turbine : str or dict
                    Name of a turbine
            hub_height : num
                    Extrapolation height

        Can also specify all of the general conversion arguments
        documented in the `convert_cutout` function.
            e.g. var_height='lml'

    """

    if "turbine" in params:
        turbine = params.pop("turbine")
        if isinstance(turbine, str):
            turbine = get_windturbineconfig(turbine)
        else:
            raise ValueError(f"Turbine ({turbine}) not found.")
        hub_height = itemgetter("hub_height")(turbine)
    elif "hub_height" in params:
        hub_height = params.pop("hub_height")
    elif "to_height" in params:
        hub_height = params.pop("to_height")
    else:
        raise ValueError("Either a turbine or hub_height must be specified.")

    params["hub_height"] = hub_height

    return cutout._convert_cutout(convert_func=convert_windspd, **params)


def windwpd(cutout: Cutout, **params):
    """
    Generate wind power density time-series

    convert.convert_cutout → convert.convert_windwpd

    Parameters
    ----------
    **params
            Must have 1 of:
                    turbine : str or dict
                            Name of a turbine
                    hub_height : num
                            Extrapolation height

            Can also specify all of the general conversion arguments
            documented in the `convert_cutout` function.
                    e.g. var_height='lml'

    """

    if "turbine" in params:
        turbine = params.pop("turbine")
        if isinstance(turbine, str):
            turbine = get_windturbineconfig(turbine)
        else:
            raise ValueError(f"Turbine ({turbine}) not found.")
        hub_height = itemgetter("hub_height")(turbine)
    elif "hub_height" in params:
        hub_height = params.pop("hub_height")
    elif "to_height" in params:
        hub_height = params.pop("to_height")
    else:
        raise ValueError("Either a turbine or hub_height must be specified.")

    params["hub_height"] = hub_height

    return cutout._convert_cutout(convert_func=convert_windwpd, **params)


def pv(
    cutout: Cutout,
    panel: str | dict,
    orientation: str | dict | Callable,
    clearsky_model: str | None = None,
    **params,
):
    """Convert downward-shortwave, upward-shortwave radiation flux and
    ambient temperature into a pv generation time-series.

    Args:
        panel (Union[str, dict]): Panel name known to the reatlas client or a panel config
            dictionary with the parameters for the electrical model in [3].
        orientation (Union[str, dict, callback]): Panel orientation can be chosen from either
            'latitude_optimal', a constant orientation {'slope': 0.0,
            'azimuth': 0.0} or a callback function with the same signature
            as the callbacks generated by the
            `geodata.pv.orientation.make_*` functions.
        clearsky_model (Optional[str]): Either the 'simple' or the 'enhanced' Reindl clearsky
            model. The default choice of None will choose dependending on
            data availability, since the 'enhanced' model also
            incorporates ambient air temperature and relative humidity.

    Returns:
        xr.DataArray: Time-series or capacity factors based on additional general
        conversion arguments.

    Note:
        You can also specify all of the general conversion arguments
        documented in the `convert_cutout` function.

    References:
        [1] Soteris A. Kalogirou. Solar Energy Engineering: Processes and Systems,
        pages 49-117,469-516. Academic Press, 2009. ISBN 0123745012.
        [2] D.T. Reindl, W.A. Beckman, and J.A. Duffie. Diffuse fraction correla-
        tions. Solar Energy, 45(1):1 - 7, 1990.
        [3] Hans Georg Beyer, Gerd Heilscher and Stefan Bofinger. A Robust Model
        for the MPP Performance of Different Types of PV-Modules Applied for
        the Performance Check of Grid Connected Systems, Freiburg, June 2004.
        Eurosun (ISES Europe Solar Congress).
    """

    if isinstance(panel, str):
        panel = get_solarpanelconfig(panel)
    if not callable(orientation):
        orientation = get_orientation(orientation)

    return cutout._convert_cutout(
        convert_func=convert_pv,
        panel=panel,
        orientation=orientation,
        clearsky_model=clearsky_model,
        **params,
    )


def pm25(cutout: Cutout, **params):
    """
    Generate PM2.5 time series 	[ug / m3]
    (see convert_pm25 for details)

    Returns:
        xr.DataArray: PM2.5 time series

    """

    return cutout._convert_cutout(convert_func=convert_pm25, **params)


__all__ = [
    "heat_demand",
    "temperature",
    "soil_temperature",
    "solar_thermal",
    "wind",
    "windspd",
    "windwpd",
    "pv",
    "pm25",
]
