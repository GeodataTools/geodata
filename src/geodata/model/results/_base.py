# Copyright 2025 Michael Davidson (UCSD), Xiqiang Liu (UCSD)

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


from __future__ import annotations

import abc
import json
import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING

import xarray as xr

from geodata.config import model_dir

if TYPE_CHECKING:
    from typing import Self

    from .._base import BaseModel

logger = logging.getLogger(__name__)


@dataclass
class BaseModelResult(abc.ABC):
    """Model result class. This class is used to store the result of a model.
    All results, regardless of the actual storage format, are grouped in monthly
    instances.

    Args:
        year (int): Year of the model.
        month (int): Month of the model.
        model (BaseModel): Model object that represent the downloaded raw data.
    """

    year: int
    month: int
    model: BaseModel

    _hashes: dict[str, str] = field(default_factory=dict)
    _prepared: bool = False
    _frequency: str = field(init=False)

    def __post_init__(self):
        if hasattr(self, "_frequency") and self._frequency != self.model.frequency:
            raise ValueError(
                f"Model frequency {self.model.frequency} does not match {self._frequency}."
            )
        self._frequency = self.model.frequency

    @property
    def frequency(self) -> str:
        """Frequency of the model."""
        return self.model.frequency

    @property
    def module(self) -> str:
        """Module of the model."""
        return self.model.source.module

    @property
    def ref_path(self) -> Path:
        """Path to the reference dataset."""
        return (
            self.model._ref_path
            / self.model.source.weather_config
            / f"{self.year:04d}"
            / f"{self.month:02d}"
        )

    @property
    def files(self) -> list[Path]:
        """List of files in the model dataset. This could potentially include any
        files that are not prepared yet."""

        atomic_files = self.model.source.get_monthly_catalog(self.year, self.month)
        if atomic_files is None:
            raise ValueError(
                f"Model files for {self.year:04d}-{self.month:02d} not found."
            )

        files = []
        for f in atomic_files:
            p = f.path.relative_to(self.ref_path)
            files.append(self.path / p)
        return files

    @property
    def ref_files(self) -> list[Path]:
        """List of files in the reference dataset."""

        atomic_files = self.model.source.get_monthly_catalog(self.year, self.month)
        if atomic_files is None:
            raise ValueError(
                f"Reference files for {self.year:04d}-{self.month:02d} not found."
            )

        return [f.path for f in atomic_files if f.path.exists()]

    @property
    def ref_params(self) -> dict[Path, Path]:
        """Reference parameters of the model."""

        ref_params = {}
        for file in self.ref_files:
            if file.name.endswith(".params.nc"):
                ref_params[file] = self.path / file.name
            else:
                ref_params[file] = self.path / f"{file.stem}.params.nc"
        return ref_params

    @property
    def path(self) -> Path:
        """Path to the model dataset."""
        return (
            model_dir
            / self.module
            / self.model.__class__.__name__
            / f"{self.year:04d}"
            / f"{self.month:02d}"
        )

    @property
    def prepared(self) -> bool:
        """Check if the model is prepared.

        Returns:
            bool: True if prepared.
        """
        if not self._prepared:
            self._prepared = self._check_prepared()
        return self._prepared

    @abc.abstractmethod
    def _check_prepared(self) -> bool:
        """Check if the model is prepared.

        Returns:
            bool: True if prepared.
        """

    @abc.abstractmethod
    def register(self, dataset: xr.Dataset):
        """Register the model result with the dataset.

        Args:
            dataset (xr.Dataset): Dataset to register.
        """

    @classmethod
    def from_year_month(cls, model: BaseModel, year: int, month: int) -> Self:
        """Create an ModelResult from a year and month.

        Args:
            model (BaseModel): BaseModel object.
            year (int): Year of the model.
            month (int): Month of the model.

        Returns:
            AtomicModel: AtomicModel object.
        """

        if not (1 <= month <= 12):
            raise ValueError(f"Month {month} is not valid. Must be between 1 and 12.")
        if not (2000 <= year <= 2100):
            raise ValueError(
                f"Year {year} is not valid. Must be between 2000 and 2100."
            )

        path = (
            model_dir
            / model.source.module
            / model.__class__.__name__
            / f"{year:04d}"
            / f"{month:02d}"
            / "meta.json"
        )

        if not path.exists():
            return cls(model=model, year=year, month=month)

        with open(path, "r") as f:
            data = json.load(f)

        if data["year"] != year or data["month"] != month:
            raise ValueError(
                f"Model year {data['year']} and month {data['month']} do not match {year} and {month}."
            )

        return cls.from_dict(data, model)

    def __repr__(self):
        return f"AtomicModel(year={self.year}, month={self.month}, ref_path={self.ref_path}, path={self.path} {len(self.files)} / {len(self.ref_files)})"

    @classmethod
    def from_dict(cls, data: dict, model: BaseModel) -> Self:
        """Create an AtomicModel from a dictionary.

        Args:
            data (dict): Dictionary with the model data.

        Returns:
            AtomicModel: AtomicModel object.
        """

        inst = cls(year=data["year"], month=data["month"], model=model)

        inst._hashes = data.get("hashes", {})
        inst.prepared
        return inst

    def to_dict(self) -> dict:
        """Convert the AtomicModel to a dictionary.

        Returns:
            dict: Dictionary with the model data.
        """

        return {
            "year": self.year,
            "month": self.month,
            "ref_path": str(self.ref_path),
            "path": str(self.path),
            "hashes": self._hashes,
            "prepared": self.prepared,
        }

    def dump(self):
        """Dump the model result to a file.

        Returns:
            dict: Dictionary with the model data.
        """
        info = self.to_dict()
        with open(self.path / "meta.json", "w") as f:
            json.dump(info, f, indent=4)

        # We dump again to check for meta preparedness
        info = self.to_dict()
        with open(self.path / "meta.json", "w") as f:
            json.dump(info, f, indent=4)

        logger.info("Model result dumped to %s", self.path / "meta.json")
