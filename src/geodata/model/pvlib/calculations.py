import numpy as np
import pandas as pd
import xarray as xr

from pvlib.atmosphere import gueymard94_pw
from pvlib.solarposition import get_solarposition



def calculate_pvlib_solarposition(ds: xr.Dataset) -> pd.DataFrame:
    """
    Wrapper for :code:`pvlib.solarposition.get_solarposition()`.  
    Allows for vectorized calculation of solar position across an xarray dataset.
    The solar zenith angle is a required input for :code:`_calculate_ghi()`.

    For full documentation on how :code:`pvlib.solarposition.get_solarposition()` calculates precipitable water,
    see: `the pvlib API reference for pvlib.solarposition.get_solarposition() <https://pvlib-python.readthedocs.io/en/v0.4.2/generated/pvlib.solarposition.get_solarposition.html>`.

    Parameters
    ----------
    ds : xarray dataset
        An xarray dataset containing series for both influx diffuse (dhi) and influx direct (dni).
    zenith : numeric
        Zenith angle of the sun in degrees, as calculated by :code:`_calculate_pvlib_solarposition()`.

    Returns
    -------
    solarposition : dataframe
        Dataframe containing solar zenith angle for a given time and set of coordinates.
    """
    nt, ny, nx = ds.sizes['time'], ds.sizes['y'], ds.sizes['x']
    time_expanded = np.broadcast_to(ds.time.values[:, None, None], (nt, ny, nx)).ravel()
    yy, xx = np.meshgrid(ds.y, ds.x, indexing="ij")
    x_expanded = np.tile(xx.ravel(), nt)
    y_expanded = np.tile(yy.ravel(), nt)
    solarposition = get_solarposition(time_expanded, y_expanded, x_expanded) # might return a pandas DataFrame or a NDArray
    multi_index = pd.MultiIndex.from_arrays([time_expanded, y_expanded, x_expanded], names=['time', 'y', 'x'])
    if isinstance(solarposition, pd.DataFrame):
        solarposition = solarposition.set_index(multi_index)
    else:
        solarposition = pd.DataFrame(solarposition)
        solarposition.index = multi_index
    return solarposition

def calculate_ghi(
    ds: xr.Dataset, 
    zenith: pd.Series
) -> xr.DataArray:
    """
    Calculates global horizontal irradiance (ghi) from data arrays representing influx diffuse (dhi) and influx direct (dni)
    Negative values are clipped.  Calculated using the formula:

    .. math::

    GHI = DHI + DNI * cos(Z)

    where Z representst the solar zenith as calculated by :code:`calculate_pvlib_solarposition()`.

    Parameters
    ----------
    ds : xarray dataset
        An xarray dataset containing series for both influx diffuse (dhi) and influx direct (dni).
    zenith : numeric
        Zenith angle of the sun in degrees, as calculated by :code:`_calculate_pvlib_solarposition()`.

    Returns
    -------
    ghi : numeric
        Global horizontal irradiance (ghi) [W m**-2].

    """
    # Convert zenith to numpy array if it's a pandas Series or xarray DataArray
    if isinstance(zenith, (pd.Series, xr.DataArray)):
        zenith_vals = zenith.values  # type: ignore[union-attr]
    else:
        zenith_vals = zenith
    # Ensure zenith_vals is a numpy array for type checking
    zenith_vals = np.asarray(zenith_vals)

    dhi = ds.influx_diffuse.values.ravel()
    dni = ds.influx_direct.values.ravel()

    if zenith_vals.size == 0:
        x_coord = ds.coords.get("x")
        if x_coord is None:
            x_coord = ds.coords.get("longitude") if "longitude" in ds.coords else ds.coords.get("lon")
        y_coord = ds.coords.get("y")
        if y_coord is None:
            y_coord = ds.coords.get("latitude") if "latitude" in ds.coords else ds.coords.get("lat")
        x_vals = np.asarray(x_coord.values) if x_coord is not None else np.array([])
        y_vals = np.asarray(y_coord.values) if y_coord is not None else np.array([])
        time_size = ds.sizes.get("time", 0)

        def _fmt_coord(arr: np.ndarray, max_show: int = 20) -> str:
            if len(arr) == 0:
                return "[]"
            if len(arr) <= max_show:
                return str(arr.tolist())
            return f"[{arr.min():g}..{arr.max():g}] (length={len(arr)})"

        raise ValueError(
            "Cannot calculate GHI: dataset has no data points. "
            "This typically occurs when xs/ys slices do not overlap with the dataset's "
            "coordinates, or when the loaded dataset is empty. "
            f"Dataset dimensions: time={time_size}, x={_fmt_coord(x_vals)}, "
            f"y={_fmt_coord(y_vals)}. "
            "Check that your xs and ys values (in degrees) overlap with these coordinate ranges."
        )

    # TODO: check if zenith is in degrees or radians and convert to radians if needed
    # it is processed from get_solarposition()
    # if zenith is in degrees, convert to radians
    if np.max(zenith_vals) > np.pi * 2:
        zenith_vals = np.deg2rad(zenith_vals)

    ghi = np.clip(
        dhi + dni * np.cos(zenith_vals),
        0,
        np.inf # `np.Inf` was removed in the NumPy 2.0 release.
    )

    reshaped_ghi = ghi.reshape(
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

    ghi.name = "ghi"
    ghi.attrs["units"] = "W m**-2"
    ghi.attrs["description"] = "Ghi calculated from influx diffuse (dhi) and influx direct (dni)."
    return ghi

def calculate_relative_humidity(
    temperature: xr.DataArray, 
    dewpoint_temperature: xr.DataArray
) -> xr.DataArray:
        """
        Calculates relative humidity based on air temperature and dewpoint temperature.
        Needed in order to calculate precipitable water using pvlib's :code:`gueymard94_pw()` function.

        Relative humidity is calculated using a version of the 
        August-Roche-Magnus equation as follows: 
        
        .. math::

            RH = 100 \cdot \frac{{\exp\left(\frac{{17.625 \cdot TD}}{{243.04 + TD}}\right)}}{{\exp\left(\frac{{17.625 \cdot T}}{{243.04 + T}}\right)}}

        where, RH is % relative humidity, TD is dew-point temperature (celsius), and T is air temperature (celsius).[#1]_ [#2]_

        Parameters
        ----------
        temperature : numeric
            Ambient air temperature at the surface. [C]
        dewpoint_temperature : numeric
            Dewpoint temperature at the surface. [C]

        Returns
        -------
        relative_humidity : numeric
            Percent relative humidity. [%]

        References
        ----------
        .. [#1] `United States Environmental Protection Agency. Hydrologic Micro Services. Meteorology - Humidity.  <https://qed.epa.gov/hms/meteorology/humidity/algorithms/>`_

        .. [#2] `University of Miami. Calculate Temperature, Dewpoint, or Relative Humidity. <https://bmcnoldy.earth.miami.edu/Humidity.html>` 

        """
        relative_humidity = 100 * (
            np.exp((17.625 * dewpoint_temperature) / (243.04 + dewpoint_temperature)) /
            np.exp((17.625 * temperature) / (243.04 + temperature))
        )
        
        # Ensure result is xarray DataArray (arithmetic operations preserve xarray types)
        if not isinstance(relative_humidity, xr.DataArray):
            relative_humidity = xr.DataArray(relative_humidity)

        relative_humidity.name = "relative_humidity"
        relative_humidity.attrs["units"] = "%"
        relative_humidity.attrs["description"] = "Relative humidity, calculated using temperature and dewpoint temperature."

        return relative_humidity

def calculate_precipitable_water(
    temperature: xr.DataArray, 
    relative_humidity: xr.DataArray
) -> xr.DataArray:
    """
    Calculates precipitable water (cm) from ambient air temperature (C) and relative humidity (%) using 
    :code:`pvlib.atmosphere.gueymard94_pw()`.  

    Precipitable water (cm) is a required input for models using CEC modules from :code:`pvlib`.
    For full documentation on how :code:`pvlib.atmosphere.gueymard94_pw()` calculates precipitable water,
    see: `the pvlib API reference for pvlib.atmosphere.gueymard94_pw() <https://pvlib-python.readthedocs.io/en/v0.4.2/generated/pvlib.atmosphere.gueymard94_pw.html>`.

    Parameters
    ----------
    temperature : numeric
        Ambient air temperature at the surface. [C]
    relative_humidity : numeric
        Percent relative humidity. [%]

    Returns
    -------
    precipitable_water : numeric
        Precipitable water (cm) calculated from ambient air temperature (C) and relative humidity (%). [cm]

    """
    # Use xarray's apply_ufunc to preserve DataArray type when calling external function
    precipitable_water = xr.apply_ufunc(
        gueymard94_pw,
        temperature,
        relative_humidity,
        dask="allowed",
        output_dtypes=[float]
    )
    precipitable_water.name = "precipitable_water"
    precipitable_water.attrs["units"] = "cm"
    precipitable_water.attrs["description"] = "Precipitable water (cm) calculated from ambient air temperature (C) and relative humidity (%)."

    return precipitable_water

def convert_kelvin_to_celsius(
    ds: xr.DataArray
) -> xr.DataArray:
    """
    Converts a temperature in Kelvin to a temperature in Celsius.

    Parameters
    ----------
    ds : numeric or xarray.DataArray
        A temperature in Kelvin [K].

    Returns
    -------
    temperature : numeric or xarray.DataArray
        A temperature in Celsius [C].
    """
    return ds - 273.15