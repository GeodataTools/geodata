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

import logging
from pathlib import Path
from typing import Optional, Union

import numpy as np
import pandas as pd
import xarray as xr
from numba import njit, prange
from tqdm.auto import tqdm

from ..config import DATASET_ROOT_PATH
from ..cutout import Cutout
from ..dataset import Dataset
from .base import BaseModel

logger = logging.getLogger(__name__)


HEIGHTS = {"u50m": 50, "u10m": 10, "u2m": 2}
__all__ = ["WindModel", "wind_speed_merra2_ext"]


@njit(parallel=True)
def _compute_wind_speed_ext(
    heights: np.ndarray, speeds: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Compute wind speed from heights and speeds.

    Args:
        heights (np.ndarray): Array of heights.
        speeds (np.ndarray): Array of speeds.

    Returns:
        tuple[np.ndarray, np.ndarray]: Tuple of coefficients (alpha, beta) and residuals.
    """
    # pylint: disable=not-an-iterable
    coeffs = np.empty(shape=heights.shape[:-1] + (2,))
    residuals = np.empty_like(speeds)
    for time in prange(heights.shape[0]):
        for lat in prange(heights.shape[1]):
            for lon in prange(heights.shape[2]):
                a = heights[time, lat, lon]
                a = np.stack((a, np.ones_like(a)), axis=-1)
                b = np.abs(speeds[time, lat, lon])

                coeff, residual, _, _ = np.linalg.lstsq(a, b)

                coeffs[time, lat, lon] = coeff
                residuals[time, lat, lon] = residual
    # pylint: enable=not-an-iterable

    return coeffs, residuals


def wind_speed_merra2_ext(
    ds: xr.Dataset,
    compute_lml: bool = True,
    half_precision: bool = True,
    compute_residuals: bool = True,
) -> xr.Dataset:
    """Compute wind speed from for MERRA2 dataset.

    Args:
        ds (xr.Dataset): MERRA2 dataset.
        compute_lml (bool, optional): Include speed at LML in calculation. Defaults to True.
        half_precision (bool, optional): Use float32 precision to store coefficients and residuals. Defaults to True.
        compute_residuals (bool, optional): Compute and store residuals against existing data points. Defaults to True.

    Returns:
        xr.Dataset: Dataset with wind speed.
    """

    disph = ds["disph"].values

    variables = [f for f in HEIGHTS if f in ds and f.replace("u", "v") in ds]
    heights = np.array([HEIGHTS[f] for f in variables]) - disph[..., np.newaxis]

    if compute_lml:
        hlml = ds["hlml"].values
        variables.append("ulml")
        heights = np.concatenate((heights, (hlml - disph)[..., np.newaxis]), axis=-1)
    log_heights = np.log(heights, out=-np.ones_like(heights), where=heights > 0)

    speeds = []
    for var in variables:
        speeds.append(((ds[var] ** 2 + ds[var.replace("u", "v")] ** 2) ** 0.5).values)
    speeds = np.stack(speeds, axis=-1).astype(heights.dtype)

    coeffs, residuals = _compute_wind_speed_ext(log_heights, speeds)

    if half_precision:
        coeffs = coeffs.astype("float32")
        residuals = residuals.astype("float32")
    ds = ds.assign_coords(coeff=["alpha", "beta"])
    ds["coeffs"] = (("time", "lat", "lon", "coeff"), coeffs)

    if compute_residuals:
        ds = ds.assign_coords(residual=variables)
        ds["residuals"] = (("time", "lat", "lon", "residual"), residuals)

    return ds[["coeffs", "residuals"]] if compute_residuals else ds[["coeffs"]]


class WindModel(BaseModel):
    """Wind model.

    Args:
        source (Union[Dataset, Cutout]): Source dataset or cutout.
        **kwargs: Keyword arguments.
    """

    type: str = "wind"

    def __init__(self, source: Union[Dataset, Cutout], **kwargs):
        super().__init__(source, **kwargs)

    def _extract_dataset_metadata(self, dataset: Dataset) -> dict:
        if not dataset.prepared:
            raise ValueError("The source dataset for this model is not prepared.")
        logger.info("Using dataset %s", dataset.module)

        metadata = {}

        metadata["name"] = metadata["module"] = dataset.module
        metadata["is_dataset"] = True
        metadata["years"] = dataset.years
        metadata["months"] = dataset.months

        return metadata

    def _extract_cutout_metadata(self, cutout: Cutout) -> dict:
        logger.info("Using cutout %s", cutout.name)

        metadata = {}

        metadata["name"] = cutout.name
        metadata["module"] = cutout.dataset_module
        metadata["years"] = cutout.years
        metadata["months"] = cutout.months

        return metadata

    def _prepare_dataset(self) -> list[tuple[str, Path]]:
        """Prepare the model from a dataset."""
        logger.info("Preparing the model from dataset.")

        prepared_files = []
        for config, file_path in tqdm(self.source.downloadedFiles):
            orig_ds_path = Path(file_path)

            ds = xr.open_dataset(orig_ds_path)
            try:
                ds = wind_speed_merra2_ext(ds)
            except SystemError:
                logger.warning(
                    "Could not compute wind speed of %s, possibly due to corrupt file.",
                    orig_ds_path.name,
                )
                continue

            ds_path = orig_ds_path.relative_to(self._ref_path).with_suffix(
                ".coeffs.nc4"
            )
            ds_path = self._path / "nc4" / ds_path
            ds_path.parent.mkdir(parents=True, exist_ok=True)
            ds.to_netcdf(ds_path)

            # (config, parameter_filepath, original_dataset_filepath)
            prepared_files.append(
                (
                    config,
                    str(ds_path.relative_to(self._path)),
                )
            )

        return prepared_files

    def _prepare_cutout(self) -> list[tuple[str, Path]]:
        """Prepare the model from a cutout."""
        logger.info("Preparing the model from dataset.")

        prepared_files = []
        for config, file_path in tqdm(self.source.downloadedFiles):
            orig_ds_path = Path(file_path)

            ds = xr.open_dataset(orig_ds_path)
            try:
                ds = wind_speed_merra2_ext(ds)
            except SystemError:
                logger.warning(
                    "Could not compute wind speed of %s, possibly due to corrupt file.",
                    orig_ds_path.name,
                )
                continue

            ds_path = orig_ds_path.relative_to(self._ref_path).with_suffix(
                ".coeffs.nc4"
            )
            ds_path = self._path / "nc4" / ds_path
            ds.to_netcdf(ds_path)

            # (config, parameter_filepath, original_dataset_filepath)
            prepared_files.append(
                (
                    config,
                    str(ds_path.relative_to(self._path)),
                )
            )

        return prepared_files

    def _estimate_dataset(
        self,
        height: int,
        years: slice,
        months: Optional[slice] = None,
        xs: Optional[slice] = None,
        ys: Optional[slice] = None,
        use_real_data: Optional[bool] = False,
    ) -> xr.Dataset:
        assert height > 0, "Height must be greater than 0."

        if months is None:
            months = slice(1, 12)

        raw_paths = [Path(p) for c, p in self.metadata["files"]]
        ds_paths = []

        for p in raw_paths:
            assert (self._path / p).exists(), f"File {p} does not exist."

            try:
                year = int(p.parts[1])
                month = int(p.parts[2])
            except ValueError:
                logger.warning(
                    "Illegal file %s contained in model, final estimated result may be incomplete",
                    p,
                )

            if year not in range(years.start, years.stop + 1) or month not in range(
                months.start, months.stop + 1
            ):
                continue

            orig_ds_path = DATASET_ROOT_PATH / self.module / p.relative_to(p.parts[0])
            ds_paths += [
                self._path / p,
                orig_ds_path.with_name(orig_ds_path.name.replace(".coeffs", "")),
            ]

        start_time = pd.Timestamp(year=years.start, month=months.start, day=1)
        end_time = pd.Timestamp(
            year=years.stop, month=months.stop + 1, day=1
        )  # exclude last month

        ds = xr.open_mfdataset(ds_paths)

        if xs is None:
            xs = ds.coords["lon"]
        if ys is None:
            ys = ds.coords["lat"]

        ds = ds.sel(lon=xs, lat=ys, time=slice(start_time, end_time))

        if height in HEIGHTS.values() and use_real_data:
            logger.info("Using real data for estimation at height %d", height)
            return (ds[f"u{height}m"] ** 2 + ds[f"v{height}m"] ** 2) ** 0.5

        alpha = ds["coeffs"][..., 0]
        beta = ds["coeffs"][..., 1]

        return alpha * np.log((height - ds["disph"]) / np.exp(-beta / alpha))

    def _estimate_cutout(self, xs: slice, ys: slice, ts: slice) -> xr.Dataset:
        raise NotImplementedError
