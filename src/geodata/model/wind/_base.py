# Copyright 2023 Michael Davidson (UCSD), Xiqiang Liu (UCSD)

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

from pathlib import Path
from typing import Callable

import xarray as xr
from tqdm.auto import tqdm

from ...logging import logger
from .._base import BaseModel

HEIGHTS = {"u50m": 50, "u10m": 10, "u2m": 2}


class WindBaseModel(BaseModel):
    """Base class for wind speed estimation (interpolation/extrapolation).
    Not meant to be instantiated directly.

    Args:
        source (Union[Dataset, Cutout]): Source dataset or cutout.
        **kwargs: Keyword arguments.
    """

    type: str = "wind"
    _prepare_fn: Callable[[xr.Dataset], xr.Dataset]

    def _prepare_dataset(self) -> list[tuple[str, Path]]:
        """Prepare the model from a dataset."""

        logger.info("Preparing the model from dataset.")

        prepared_files = []
        for file_path in tqdm(self.metadata["files_orig"], dynamic_ncols=True):
            orig_ds_path: Path = self._ref_path / file_path
            ds = xr.open_dataset(orig_ds_path)
            try:
                ds = self._prepare_fn(ds)
            except SystemError:
                logger.warning(
                    "Could not compute wind speed of %s, possibly due to corrupt file.",
                    orig_ds_path.name,
                )
                continue

            ds_path: Path = (
                self._path / "nc4" / Path(file_path).with_suffix(".params.nc4")
            )
            ds_path.parent.mkdir(parents=True, exist_ok=True)
            ds.to_netcdf(ds_path)

            prepared_files.append(str(ds_path.relative_to(self._path)))

        return prepared_files

    def _prepare_cutout(self) -> list[tuple[str, Path]]:
        """Prepare the model from a cutout."""

        logger.info("Preparing the model from cutout.")
        prepared_files = []

        for yearmonth in tqdm(self.source.coords["year-month"].to_index()):
            orig_ds_path = Path(self.source.datasetfn(yearmonth))

            ds = xr.open_dataset(orig_ds_path)
            try:
                ds = self._prepare_fn(ds)
            except SystemError:
                logger.warning(
                    "Could not compute wind speed of %s, possibly due to corrupt file.",
                    orig_ds_path.name,
                )
                continue

            ds_path = orig_ds_path.relative_to(self._ref_path).with_suffix(
                ".params.nc4"
            )
            ds_path = self._path / "nc4" / ds_path
            ds_path.parent.mkdir(parents=True, exist_ok=True)
            ds.to_netcdf(ds_path)

            prepared_files.append(str(ds_path.relative_to(self._path)))

        return prepared_files
