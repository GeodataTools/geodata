# Copyright 2023, 2025 Michael Davidson (UCSD), Xiqiang Liu (UCSD)

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

"""Wind speed estimation based on the Model interface.

This module provides a wind speed estimation model based on the Model interface.

Extrapolation:
    The wind speed is estimated using a linear regression of the logarithm of the wind speed
    against the logarithm of the height. The coefficients of the regression are stored in the
    dataset as a 2D array of shape (2,). The first coefficient is the slope of the regression
    line, and the second coefficient is the intercept. The wind speed is then estimated as
    :math:`\\hat{v} = \\alpha \\cdot \\log(z) + \\beta`, where :math:`\\alpha` and :math:`\\beta`
    are the coefficients, and :math:`z` is the height.

Residuals:
    The residuals of the regression are stored in the dataset as a 2D array of shape (n,).
    The residuals are the difference between the estimated wind speed and the actual wind speed.
    The residuals are stored in the same order as the heights used in the regression.

Example:
    >>> from geodata import Cutout
    >>> from geodata.model import WindModel
    >>> cutout = Cutout("merra2", 2010, 2010, 0, 0, 0, 0)
    >>> model = WindModel(cutout)
    >>> model.prepare()
    >>> model.estimate(xs=slice(1, 2), ys=slice(1, 2), years=slice(2010, 2010), months=slice(1, 2))
"""

import xarray as xr

from ...resource import get_windturbineconfig
from .._base import BaseModel

from scipy.interpolate import interp1d

HEIGHTS = {"u50m": 50, "u10m": 10, "u2m": 2}


class WindBaseModel(BaseModel):
    """Base class for wind speed estimation (interpolation/extrapolation).
    Not meant to be instantiated directly.

    Args:
        source (Union[Dataset, Cutout]): Source dataset or cutout.
        **kwargs: Keyword arguments.
    """

    type: str = "wind"

    def estimate(
        self,
        years: slice | None = None,
        months: slice | None = None,
        xs: slice | None = None,
        ys: slice | None = None,
        **kwargs,
    ) -> xr.DataArray:
        """Estimate wind speed or CF at the given locations and times. If a turbine is
        specified, the CF is calculated based on the wind speed and the turbine's power curve.
        Otherwise, pass in the `height` keyword argument to estimate wind speed at a specific height.

        Args:
            years (slice, optional): Years.
            months (slice, optional): Months. If None, all months are estimated.
            xs (slice, optional): X coordinates. If None, all x coordinates in source are estimated.
            ys (slice, optional): Y coordinates. If None, all y coordinates in source are estimated.
            **kwargs: Additional keyword arguments to pass to the model.
                - `turbine` (str): Name of the wind turbine to estimate power output.
                - `height` (int): Height at which to estimate wind speed. If not specified, the model will use the default height.

        Raises:
            ValueError: If neither 'turbine' nor 'height' is specified in kwargs.

        Returns:
            xr.DataArray: Estimated wind speed.
        """

        if "turbine" in kwargs:
            return self._estimate_power(
                years=years, months=months, xs=xs, ys=ys, **kwargs
            )

        if "height" not in kwargs:
            raise ValueError(
                "Either 'turbine' or 'height' must be specified to estimate wind speed."
            )

        return super().estimate(years, months, xs, ys, **kwargs)

    def _estimate_power(
        self,
        turbine: str,
        xs: slice | None = None,
        ys: slice | None = None,
        years: slice | None = None,
        months: slice | None = None,
    ) -> None:
        """Estimate wind speed at the given locations and times.

        Args:
            turbine (str): Turbine name.
            xs (slice, optional): X slice. Defaults to None.
            ys (slice, optional): Y slice. Defaults to None.
            years (slice, optional): Year slice. Defaults to None.
            months (slice, optional): Month slice. Defaults to None.

        Returns:
            xr.DataArray: Estimated wind speed.
        """

        # Get the wind turbine configuration
        try:
            turbineconf = get_windturbineconfig(turbine)
        except FileNotFoundError:
            raise ValueError(f"Wind turbine configuration '{turbine}' not found.")

        speed = self.estimate(
            years=years, months=months, xs=xs, ys=ys, height=turbineconf["hub_height"]
        )

        interp_fn = interp1d(
            turbineconf["V"],
            turbineconf["POW"],
            bounds_error=False,
            fill_value="extrapolate",
        )

        # Calculate the power output
        power = xr.apply_ufunc(
            interp_fn,
            speed,
            vectorize=True,
            dask="parallelized",
            output_dtypes=[float],
        )

        return xr.Dataset({"cf": power / turbineconf["P"]})
