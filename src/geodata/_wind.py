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

import numpy as np

logger = logging.getLogger(__name__)


"""
Extrapolation functions
	Note: these are called internally
"""


def log_ratio(ds, to_height, from_height, from_name):
    """Logarithmic ratio law
    # Equation (2) in Andresen, G. et al (2015):
    # 	'Validation of Danish wind time series from a new global renewable
    # 	energy atlas for energy system analysis'.
    """

    wnd_spd = ds[from_name] * (
        np.log(to_height / ds["roughness"]) / np.log(ds[from_height] / ds["roughness"])
    )
    wnd_spd.attrs.update(
        {
            "long_name": f"extrapolated {to_height} m wind speed using log ratio",
            "units": "m s**-1",
        }
    )
    return wnd_spd


def log_law(ds, to_height, from_height, from_name):
    """Logarithmic (integration) law
    [3] S. Emeis, Wind Energy Meteorology (Springer, Berlin, 2013).
    """
    vonk = 0.4
    wnd_spd = ds[from_name] + (
        ds["ustar"] / vonk * np.log((to_height - ds["disph"]) / ds[from_height])
    )
    wnd_spd.attrs.update(
        {
            "long_name": f"extrapolated {to_height} m wind speed using log (integration) law",
            "units": "m s**-1",
        }
    )
    return wnd_spd


## Flux stability correction functions
def psi_linear(z, ds):
    """NOTE: This function does not perform well for high z/L
    Linear stability correction [1,2]
            z = to_height
            ds = dataset with variables: L

    Eval = 	0, if z/L <=0
                    linear if z/L > 0
    [1] Businger, J.A., Wyngaard, J.C., Izumi, Y., Bradley, E.F., 1971. Flux profile relationships in the atmospheric surface layer. J. Atmos. Sci. 28, 181–189.
    [2] Dyer, A.J., 1974. A review of flux–profile relationships. Boundary-Layer Meteorol. 7, 363–372.
    """
    beta = 5.2
    ds["a"] = z / ds["L"]
    ds["psim"] = ds["a"] * 0  # create zeroes same length as dataset
    ds["psim"].values[0 < ds["a"]] = -beta * ds["a"].values[0 < ds["a"]]
    ds["psim"].values[ds["a"] <= 0] = 0
    return ds["psim"]


def psi_linearexp(z, ds):
    """Linear-exponential piecewise stability correction [1] (repeated in [2])
            z = to_height
            ds = dataset with variables: L

    [1] Emeis, S. (2013). Wind Energy Meteorology. Retrieved from http://link.springer.com/10.1007/978-3-642-30523-8 (Note: error in eq 3.21)
    [2] Rose, S., & Apt, J. (2016). Quantifying sources of uncertainty in reanalysis derived wind speed. Renewable Energy, 94, 157–165.
            https://doi.org/10.1016/j.renene.2016.03.028
    """
    aconst = 5
    A = 1
    B = 2 / 3
    C = 5
    D = 0.35
    ds["a"] = z / ds["L"]
    ds["psim"] = ds["a"] * 0  # create zeroes same length as dataset
    ds["psim"].values[(0 < ds["a"]) & (ds["a"] <= 0.5)] = (
        -aconst * ds["a"].values[(0 < ds["a"]) & (ds["a"] <= 0.5)]
    )
    ds["psim"].values[0.5 < ds["a"]] = -A * (
        ds["a"].values[0.5 < ds["a"]]
        + B
        * (ds["a"].values[0.5 < ds["a"]] - C / D)
        * np.exp(-D * ds["a"].values[0.5 < ds["a"]])
        + B * C / D
    )
    ds["psim"].values[ds["a"] <= 0] = 0
    return ds["psim"]


def psi_linearexpconst(z, ds, const=7):
    """Linear-exponential piecewise stability correction [1] (repeated in [2]) with constant plateau
            z = to_height
            ds = dataset with variables: L
            const = upper bound of z/L, after which = constant

    [1] Emeis, S. (2013). Wind Energy Meteorology. Retrieved from http://link.springer.com/10.1007/978-3-642-30523-8 (Note: error in eq 3.21)
    [2] Rose, S., & Apt, J. (2016). Quantifying sources of uncertainty in reanalysis derived wind speed. Renewable Energy, 94, 157–165.
            https://doi.org/10.1016/j.renene.2016.03.028
    """
    aconst = 5
    A = 1
    B = 2 / 3
    C = 5
    D = 0.35
    ds["a"] = z / ds["L"]
    ds["psim"] = ds["a"] * 0  # create zeroes same length as dataset
    ds["psim"].values[(0 < ds["a"]) & (ds["a"] <= 0.5)] = (
        -aconst * ds["a"].values[(0 < ds["a"]) & (ds["a"] <= 0.5)]
    )
    ds["psim"].values[0.5 < ds["a"]] = -A * (
        ds["a"].values[0.5 < ds["a"]]
        + B
        * (ds["a"].values[0.5 < ds["a"]] - C / D)
        * np.exp(-D * ds["a"].values[0.5 < ds["a"]])
        + B * C / D
    )
    ds["psim"].values[ds["a"] > const] = -A * (
        const + B * (const - C / D) * np.exp(-D * const) + B * C / D
    )
    ds["psim"].values[ds["a"] <= 0] = 0
    return ds["psim"]


def L_vph(ds):
    """Obuhkov length using virtual potential heat flux term [1] (described in detail in SI [2])

    [1] Emeis, S. (2013). Wind Energy Meteorology. Retrieved from http://link.springer.com/10.1007/978-3-642-30523-8 (Note: error in eq 3.21)
    [2] Rose, S., & Apt, J. (2016). Quantifying sources of uncertainty in reanalysis derived wind speed. Renewable Energy, 94, 157–165.
            https://doi.org/10.1016/j.renene.2016.03.028
    """
    vonk = 0.4  # Von Karman constant
    grav = 9.81  # gravitational acceleration in kg m s-2
    CPD = 1004  # specific heat of dry air at constant pressure J K-1 kg-1
    Le = 2.257e6  # latent heat of evaporation [J/kg]
    kp = 2 / 7  # Poisson constant
    Rd = 287  #  Ideal gas constant [J/kg/K]
    p0 = 1e5  # standard air pressure

    ds["p"] = ds["rhoa"] * Rd * ds["tlml"]
    ds["vphflux"] = (
        ds["hflux"] + 0.61 * CPD / Le * ds["tlml"] * (p0 / ds["p"]) ** kp * ds["eflux"]
    )
    ds["L"] = -(ds["tlml"] * ds["ustar"] ** 3 * CPD * ds["rhoa"]) / (
        vonk * grav * ds["vphflux"]
    )
    return ds["L"]


def winddir(ds):
    """Wind direction using lowest model layer"""

    ds["winddir"] = np.degrees(np.arctan(ds["ulml"] / ds["vlml"]))
    ds["winddir"].values[ds["vlml"] < 0] += 180
    ds["winddir"].values[(ds["vlml"] > 0) & (ds["ulml"] < 0)] += 360
    return ds["winddir"]


def _log_law_flux(ds, to_height, from_height, from_name, psifn, Lfn=L_vph):
    """Compute logarithmic (integration) law given stability correction fn in terms of Obukhov length (derived from heat flux) [1]
    Called by: log_law_flux_**

    [1] Sharan, M., & Aditi. (2009).
            Performance of various similarity functions for nondimensional wind and temperature profiles in the surface layer in stable conditions.
            Atmospheric Research, 94(2), 246-253.
            https://doi.org/10.1016/j.atmosres.2009.05.014
    """
    vonk = 0.4  # Von Karman constant
    ds["L"] = L_vph(ds)

    wnd_spd = ds[from_name] + ds["ustar"] / vonk * (
        np.log((to_height - ds["disph"]) / ds[from_height])
        - psifn(to_height, ds[["L", "roughness"]])
    )
    wnd_spd.attrs.update(
        {
            "long_name": f"extrapolated {to_height} m wind speed using log (integration) law and stability correction {psifn}",
            "units": "m s**-1",
        }
    )
    return wnd_spd


def log_law_flux_linear(ds, to_height, from_height, from_name):
    """Logarithmic (integration) law with linear stability correction in terms of Obukhov length"""
    return _log_law_flux(ds, to_height, from_height, from_name, psi_linear)


def log_law_flux_linearexp(ds, to_height, from_height, from_name):
    """Logarithmic (integration) law with piecewise linear-exponential stability correction in terms of Obukhov length"""
    return _log_law_flux(ds, to_height, from_height, from_name, psi_linearexp)


def log_law_flux_linearexpconst(ds, to_height, from_height, from_name):
    """Logarithmic (integration) law with piecewise linear-exponential-constant stability correction in terms of Obukhov length"""
    return _log_law_flux(ds, to_height, from_height, from_name, psi_linearexpconst)


"""
Main call (from convert.convert_wind)
"""


def extrapolate_wind_speed(
    ds, to_height, extrap_fn=log_ratio, from_height=None, var_height=None
):
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

    to_name = "wnd{h:0d}m".format(h=int(to_height))
    if to_name in ds:
        # already found wind speed at given height in dataset
        return ds[to_name]

    # Sanitize roughness for logarithm: 0.0002 corresponds to open water [2]
    ds["roughness"].values[ds["roughness"].values <= 0.0] = 0.0002

    if from_height is not None:
        # passed a from_height
        if var_height is not None:
            raise AssertionError(
                "Cannot pass both from_height and var_height to extrapolate_wind_speed"
            )
        from_name = "wnd{h:0d}m".format(h=int(from_height))
        ds["from_height"] = from_height

        wnd_spd = extrap_fn(ds, to_height, "from_height", from_name)
        wnd_spd.attrs["long_name"] = (
            wnd_spd.attrs["long_name"] + ", " + f"from fixed height = {from_height}"
        )

    elif var_height is not None:
        # passed a variable height (eg lml)
        # set variable names
        from_height = f"h{var_height}"
        from_name = f"wnd{var_height}"

        wnd_spd = extrap_fn(ds, to_height, from_height, from_name)
        wnd_spd.attrs["long_name"] = (
            wnd_spd.attrs["long_name"] + ", " + f"from variable height = {var_height}"
        )

    else:
        # based on nearest height
        heights = np.asarray([int(s[3:-1]) for s in ds if s.startswith("wnd")])
        if len(heights) == 0:
            raise AssertionError("Wind speed is not in dataset")

        from_height = heights[np.argmin(np.abs(heights - to_height))]
        from_name = "wnd{h:0d}m".format(h=int(from_height))
        ds["from_height"] = from_height

        wnd_spd = extrap_fn(ds, to_height, "from_height", from_name)
        wnd_spd.attrs["long_name"] = (
            wnd_spd.attrs["long_name"] + ", " + f"from nearest height = {from_height}"
        )

    return wnd_spd.rename(to_name)
