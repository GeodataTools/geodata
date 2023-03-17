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
# along with this program. If not, see <http://www.gnu.org/licenses/>.


"""
GEODATA

Geospatial Data Collection and "Pre-Analysis" Tools
"""

from __future__ import absolute_import

import datetime as dt
import logging
from operator import itemgetter

import numpy as np
import xarray as xr
from six import string_types

from . import wind as windm
from .pv.irradiation import TiltedIrradiation
from .pv.orientation import SurfaceOrientation, get_orientation
from .pv.solar_panel_model import SolarPanelModel
from .pv.solar_position import SolarPosition
from .resource import get_solarpanelconfig, get_windturbineconfig, windturbine_smooth
from .utils import make_optional_progressbar

logger = logging.getLogger(__name__)


def convert_cutout(cutout, convert_func, show_progress=False, **convert_kwds):
    """
    Convert and aggregate a weather-based renewable generation time-series.

    NOTE: Not meant to be used by the user him or herself. Rather it is a
    gateway function that is called by all the individual time-series
    generation functions like pv and wind. Thus, all its parameters are also
    available from these.

    Parameters (passed through as **params)
    ---------------------------------------
    show_progress : boolean|string
            Whether to show a progress bar if boolean and its label if given as a
            string (defaults to True).

    Returns
    -------
    resource : xr.DataArray
            Time-series of renewable generation aggregated to buses, if
            `matrix` or equivalents are provided else the total sum of
            generated energy.

    Internal Parameters (provided by f.ex. wind and pv)
    ---------------------------------------------------
    convert_func : Function
            Callback like convert_wind, convert_pv
    """
    assert cutout.prepared, "The cutout has to be prepared first."

    results = []

    yearmonths = cutout.coords["year-month"].to_index()

    if isinstance(show_progress, string_types):
        prefix = show_progress
    else:
        func_name = (
            convert_func.__name__[len("convert_") :]
            if convert_func.__name__.startswith("convert_")
            else convert_func.__name__
        )
        prefix = "Convert `{}`: ".format(func_name)

    maybe_progressbar = make_optional_progressbar(
        show_progress, prefix, len(yearmonths)
    )

    for ym in maybe_progressbar(yearmonths):
        with xr.open_dataset(cutout.datasetfn(ym)) as ds:
            if "view" in cutout.meta.attrs:
                if isinstance(cutout.meta.attrs["view"], str):
                    cutout.meta.attrs["view"] = {}
                    cutout.meta.attrs.setdefault("view", {})["x"] = slice(
                        min(cutout.meta.coords["x"]).values.tolist(),
                        max(cutout.meta.coords["x"]).values.tolist(),
                    )
                    cutout.meta.attrs.setdefault("view", {})["y"] = slice(
                        min(cutout.meta.coords["y"]).values.tolist(),
                        max(cutout.meta.coords["y"]).values.tolist(),
                    )
                ds = ds.sel(**cutout.meta.attrs["view"])

            da = convert_func(ds, **convert_kwds)
            results.append(da.load())

    results = xr.concat(results, dim="time")

    return results


## temperature


def convert_temperature(ds):
    """Return outside temperature (useful for e.g. heat pump T-dependent
    coefficient of performance).
    """

    # Temperature is in Kelvin
    return ds["temperature"] - 273.15


def temperature(cutout, **params):
    return cutout.convert_cutout(convert_func=convert_temperature, **params)


## soil temperature


def convert_soil_temperature(ds):
    """Return soil temperature (useful for e.g. heat pump T-dependent
    coefficient of performance).
    """

    # Temperature is in Kelvin

    # There are nans where there is sea; by setting them
    # to zero we guarantee they do not contribute when multiplied
    # by matrix in geodata/aggregate.py
    return (ds["soil temperature"] - 273.15).fillna(0.0)


def soil_temperature(cutout, **params):
    return cutout.convert_cutout(convert_func=convert_soil_temperature, **params)


## heat demand


def convert_heat_demand(ds, threshold, a, constant, hour_shift):
    # Temperature is in Kelvin; take daily average
    T = ds["temperature"]
    T.coords["time"].values += np.timedelta64(dt.timedelta(hours=hour_shift))

    T = ds["temperature"].resample(time="1D").mean(dim="time")
    threshold += 273.15
    heat_demand_value = a * (threshold - T)

    heat_demand_value.values[heat_demand_value.values < 0.0] = 0.0

    return constant + heat_demand_value


def heat_demand(cutout, threshold=15.0, a=1.0, constant=0.0, hour_shift=0.0, **params):
    """
    Convert outside temperature into daily heat demand using the
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

    Parameters
    ----------
    threshold : float
            Outside temperature in degrees Celsius above which there is no
            heat demand.
    a : float
            Linear factor relating heat demand to outside temperature.
    constant : float
            Constant part of heat demand that does not depend on outside
            temperature (e.g. due to water heating).
    hour_shift : float
            Time shift relative to UTC for taking daily average

    Note
    ----
    You can also specify all of the general conversion arguments
    documented in the `convert_cutout` function.
    """

    return cutout.convert_cutout(
        convert_func=convert_heat_demand,
        threshold=threshold,
        a=a,
        constant=constant,
        hour_shift=hour_shift,
        **params,
    )


## solar thermal collectors


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


def solar_thermal(
    cutout,
    orientation=None,
    trigon_model="simple",
    clearsky_model="simple",
    c0=0.8,
    c1=3.0,
    t_store=80.0,
    **params,
):
    """
    Convert downward short-wave radiation flux and outside temperature
    into time series for solar thermal collectors.

    Mathematical model and defaults for c0, c1 based on model in [1].

    Parameters
    ----------
    cutout : cutout
    orientation : dict or str or function
            Panel orientation with slope and azimuth (units of degrees), or
            'latitude_optimal'.
    trigon_model : str
            Type of trigonometry model
    clearsky_model : str or None
            Type of clearsky model for diffuse irradiation. Either
            `simple' or `enhanced'.
    c0, c1 : float
            Parameters for model in [1] (defaults to 0.8 and 3., respectively)
    t_store : float
            Store temperature in degree Celsius

    Note
    ----
    You can also specify all of the general conversion arguments
    documented in the `convert_cutout` function.

    References
    ----------
    [1] Henning and Palzer, Renewable and Sustainable Energy Reviews 30
            (2014) 1003-1018
    """
    if orientation is None:
        orientation = {"slope": 45.0, "azimuth": 180.0}

    if not callable(orientation):
        orientation = get_orientation(orientation)

    return cutout.convert_cutout(
        convert_func=convert_solar_thermal,
        orientation=orientation,
        trigon_model=trigon_model,
        clearsky_model=clearsky_model,
        c0=c0,
        c1=c1,
        t_store=t_store,
        **params,
    )


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


def wind(cutout, turbine, smooth=False, **params):
    """
    Generate wind generation time-series

    - loads turbine dict based on passed parameters  	(resource.get_windturbineconfig)
    - optionally, smooths turbine power curve 			(resource.windturbine_smooth)
    - calls convert_wind								(convert.convert_cutout)

    Parameters
    ----------
    turbine : str or dict
            Name of a turbine known by the reatlas client or a
            turbineconfig dictionary with the keys 'hub_height' for the
            hub height and 'V', 'POW' defining the power curve.
    smooth : bool or dict
            If True smooth power curve with a gaussian kernel as
            determined for the Danish wind fleet to Delta_v = 1.27 and
            sigma = 2.29. A dict allows to tune these values.

    Note
    ----
    You can also specify all of the general conversion arguments
    documented in the `convert_cutout` function.

    References
    ----------
    [1] Andresen G B, Søndergaard A A and Greiner M 2015 Energy 93, Part 1
            1074 – 1088. doi:10.1016/j.energy.2015.09.071
    """

    if isinstance(turbine, string_types):
        turbine = get_windturbineconfig(turbine)

    if smooth:
        turbine = windturbine_smooth(turbine, params=smooth)

    return cutout.convert_cutout(convert_func=convert_wind, turbine=turbine, **params)


def windspd(cutout, **params):
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
        if isinstance(turbine, string_types):
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

    return cutout.convert_cutout(convert_func=convert_windspd, **params)


def windwpd(cutout, **params):
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
        if isinstance(turbine, string_types):
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

    return cutout.convert_cutout(convert_func=convert_windwpd, **params)


## solar PV


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


def pv(cutout, panel, orientation, clearsky_model=None, **params):
    """
    Convert downward-shortwave, upward-shortwave radiation flux and
    ambient temperature into a pv generation time-series.

    Parameters
    ----------
    panel : str or dict
            Panel name known to the reatlas client or a panel config
            dictionary with the parameters for the electrical model in [3].
    orientation : str, dict or callback
            Panel orientation can be chosen from either
            'latitude_optimal', a constant orientation {'slope': 0.0,
            'azimuth': 0.0} or a callback function with the same signature
            as the callbacks generated by the
            `geodata.pv.orientation.make_*' functions.
    clearsky_model : str or None
            Either the 'simple' or the 'enhanced' Reindl clearsky
            model. The default choice of None will choose dependending on
            data availability, since the 'enhanced' model also
            incorporates ambient air temperature and relative humidity.

    Returns
    -------
    pv : xr.DataArray
            Time-series or capacity factors based on additional general
            conversion arguments.

    Note
    ----
    You can also specify all of the general conversion arguments
    documented in the `convert_cutout` function.

    References
    ----------
    [1] Soteris A. Kalogirou. Solar Energy Engineering: Processes and Systems,
            pages 49–117,469–516. Academic Press, 2009. ISBN 0123745012.
    [2] D.T. Reindl, W.A. Beckman, and J.A. Duffie. Diffuse fraction correla-
            tions. Solar Energy, 45(1):1 – 7, 1990.
    [3] Hans Georg Beyer, Gerd Heilscher and Stefan Bofinger. A Robust Model
            for the MPP Performance of Different Types of PV-Modules Applied for
            the Performance Check of Grid Connected Systems, Freiburg, June 2004.
            Eurosun (ISES Europe Solar Congress).
    """

    if isinstance(panel, string_types):
        panel = get_solarpanelconfig(panel)
    if not callable(orientation):
        orientation = get_orientation(orientation)

    return cutout.convert_cutout(
        convert_func=convert_pv,
        panel=panel,
        orientation=orientation,
        clearsky_model=clearsky_model,
        **params,
    )


## Air pollution


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
    Evaluation of the surface PM2.5 in Version 1 of the NASA MERRA Aerosol Reanalysis over the United States. Atmospheric Environment, 125, 100–111.
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


def pm25(cutout, **params):
    """
    Generate PM2.5 time series 	[ug / m3]
    (see convert_pm25 for details)

    Parameters
    ----------
    **params : None needed currently.

    Returns
    -------
    pm25 : xr.DataArray

    """

    return cutout.convert_cutout(convert_func=convert_pm25, **params)


## Manipulate arbitrary variables


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
