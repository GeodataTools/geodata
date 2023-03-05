## Copyright 2023 Michael Davidson (UCSD), Xiqiang Liu (UCSD)

## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License as
## published by the Free Software Foundation; either version 3 of the
## License, or (at your option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.

## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.


import abc
import logging
from dataclasses import dataclass
from typing import Union

import xarray as xr

import geodata

logger = logging.getLogger(__name__)


@dataclass
class _Source:
    """A unified representation of source data for a model.

    Args:
        name (str): The name of the source.
        source (Union[geodata.Dataset, geodata.Cutout, xr.Dataset]): The source of the model. Can be a Dataset, Cutout or xarray.Dataset.
        interpolate (bool, optional): Interpolate the source to the same grid as the target. Defaults to False.
        **kwargs: Additional keyword arguments to pass to the model.
    """

    name: str
    source: Union[geodata.Dataset, geodata.Cutout, xr.Dataset]
    interpolate: bool = False


class BaseModel(abc.ABC):
    """Base class for geospatial modeling.

    Args:
        name (str): The name of the model.
        source (Union[geodata.Dataset, geodata.Cutout, xr.Dataset]): The source of the model. Can be a Dataset, Cutout or xarray.Dataset.
        interpolate (bool, optional): Interpolate the source to the same grid as the target. Defaults to False.
        **kwargs: Additional keyword arguments to pass to the model.
    """

    def __init__(
        self,
        name: str,
        source: Union[geodata.Dataset, geodata.Cutout, xr.Dataset],
        interpolate: bool = False,
        **kwargs,
    ):
        self.name = name
        self.source = self._prepare(source)
        self.interpolate = interpolate
        self.kwargs = kwargs

    def __repr__(self):
        return f"Model(name={self.name}, source={self.source}, interpolate={self.interpolate})"

    @abc.abstractmethod
    def estimate(self, xs: slice, ys: slice, ts: slice) -> xr.Dataset:
        """Estimate wind speed from the source data.

        Args:
            xs (slice): The x slice to estimate.
            ys (slice): The y slice to estimate.
            ts (slice): The t slice to estimate.

        Returns:
            xr.Dataset: The estimated wind speed.
        """

    @abc.abstractmethod
    def _prepare(
        self, source: Union[geodata.Dataset, geodata.Cutout, xr.Dataset]
    ) -> _Source:
        """Prepare the source data.

        Args:
            source (Union[geodata.Dataset, geodata.Cutout, xr.Dataset]): The source of the model. Can be a Dataset, Cutout or xarray.Dataset.

        Returns:
            _Source: The prepared source.
        """
