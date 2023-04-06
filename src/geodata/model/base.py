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


import abc
import json
import logging
import pickle
import shutil
from pathlib import Path
from typing import Union

import xarray as xr

import geodata

from ..config import model_dir

logger = logging.getLogger(__name__)


class ModelMetaSerializer(json.encoder.JSONEncoder):
    r"""JSON encoder for model metadata that serializes `slice(...)` as `str`."""

    def default(self, o):
        if isinstance(o, slice):
            return str(o)
        return json.encoder.JSONEncoder.default(self, o)


class BaseModel(abc.ABC):
    """Base class for geospatial modeling.

    Args:
        name (str): The name of the model.
        source (Union[geodata.Dataset, geodata.Cutout, xr.Dataset]): The source of the model. Can be a Dataset, Cutout or xarray.Dataset.
        interpolate (bool, optional): Interpolate the source to the same grid as the target. Defaults to False.
        **kwargs: Additional keyword arguments to pass to the model.
    """

    def __init__(self, source: Union[geodata.Dataset, geodata.Cutout], **kwargs):
        if not source.prepared:
            raise ValueError(
                "The source Dataset/Cutout for this model is not prepared."
            )

        self.source = source
        self._extra_kwargs = kwargs

        if isinstance(source, geodata.Dataset):
            self._ref_path = model_dir.parent / self.source.module
            self._path = model_dir / self.type / self.source.module
        else:
            self._ref_path = model_dir.parent / "cutouts"
            self._path = model_dir / self.type / self.source.name

        if self._check_prepared():
            with open(self._path / "meta.json", encoding="utf_8") as f:
                self._metadata_json = json.load(f)
            with open(self._path / "meta.pkl", "rb") as f:
                self._metadata = pickle.load(f)
        else:
            self.metadata = source

    def __repr__(self):
        return f"Model(source={self.source}, type={self.type})"

    @property
    def name(self) -> str:
        """Name of the model."""
        return f"{self._metadata['name']}_{self.type}"

    @property
    def module(self) -> str:
        """Module of the model."""
        return self._metadata["name"]

    @property
    @abc.abstractmethod
    def type(self) -> str:
        """Type of the model."""

    def estimate(
        self, *args: Union[slice, int], **kwargs: Union[slice, int]
    ) -> xr.Dataset:
        """Estimate the wind speed at given coordinates.

        Args:
            height (int): Height of the wind speed, need to be greater than 0.
            years (slice): Years.
            months (slice, optional): Months. If None, all months are estimated.
            xs (slice): X coordinates. If None, all x coordinates in source are estimated.
            ys (slice): Y coordinates. If None, all y coordinates in source are estimated.
            use_real_data (bool, optional): If available, use real data for estimation. Defaults to False.

        Returns:
            xr.Dataset: Dataset with wind speed.
        """

        if self.from_dataset:
            return self._estimate_dataset(*args, **kwargs)
        else:
            return self._estimate_cutout(*args, **kwargs)

    @property
    def metadata(self) -> dict:
        """The metadata of the model."""
        if hasattr(self, "_metadata"):
            return self._metadata
        else:
            raise ValueError("Metadata has not be prepared yet.")

    @metadata.setter
    def metadata(self, source: Union[geodata.Dataset, geodata.Cutout]):
        if isinstance(source, geodata.Dataset):
            self._metadata = self._extract_dataset_metadata(source)
        elif isinstance(source, geodata.Cutout):
            self._metadata = self._extract_cutout_metadata(source)
        else:
            raise ValueError(
                "The model is instantiated without a valid source, such as a Dataset or Cutout."
            )

    @abc.abstractmethod
    def _extract_dataset_metadata(self, dataset: geodata.Dataset) -> dict:
        """Extract metadata from dataset.

        Args:
            dataset (Dataset): Dataset.

        Returns:
            dict: Metadata.
        """

    @abc.abstractmethod
    def _extract_cutout_metadata(self, cutout: geodata.Cutout) -> dict:
        """Extract metadata from cutout.

        Args:
            cutout (Cutout): Cutout.

        Returns:
            dict: Metadata.
        """

    @property
    def prepared(self) -> bool:
        """Check if the model is prepared.

        Returns:
            bool: True if prepared.
        """
        return self._check_prepared()

    def _check_prepared(self) -> bool:
        """Check if the model is prepared.

        Returns:
            bool: True if prepared.
        """

        assert self._path is not None, "The model saving path has not been set yet."

        nc4_path: Path = self._path / "nc4"
        meta_path: Path = self._path / "meta"

        return (
            nc4_path.exists()
            and meta_path.with_suffix(".pkl").exists()
            and meta_path.with_suffix(".json").exists()
            and len(list(nc4_path.rglob("*.nc4"))) == len(self.source.downloadedFiles)
        )

    def prepare(self, force: bool = False):
        """Prepare the model.

        Args:
            force (bool, optional): Force re-prepare the model. Defaults to False."""

        if self.prepared and not force:
            logger.info("The model is already prepared.")
            return

        if not self._check_prepared():
            logger.info("Model not present in model directory, creating.")
            shutil.rmtree(self._path, ignore_errors=True)
            (self._path / "nc4").mkdir(exist_ok=True, parents=True)

        if self.from_dataset:
            self._metadata["files"] = self._prepare_dataset()
        else:
            self._metadata["files"] = self._prepare_cutout()

        with open(self._path / "meta.json", "w", encoding="utf-8") as f:
            json.dump(self._metadata, f, cls=ModelMetaSerializer, indent=4)

        with open(self._path / "meta.pkl", "wb") as f:
            pickle.dump(self._metadata, f)

        logger.info("Finished preparing model.")

    @abc.abstractmethod
    def _prepare_dataset(self) -> list:
        """Prepare the model from a dataset.

        Returns:
            list: List of files.
        """

    @abc.abstractmethod
    def _prepare_cutout(self) -> list:
        """Prepare the model from a cutout.

        Returns:
            list: List of files.
        """

    @property
    def from_dataset(self) -> bool:
        """Check if the model is from a dataset."""
        return self._metadata["is_dataset"]

    @abc.abstractmethod
    def _estimate_dataset(self, *args, **kwargs) -> xr.Dataset:
        """Estimate the wind speed from a dataset.

        Args:
            height (int): Height of the wind speed, need to be greater than 0.
            years (slice): Years.
            months (slice, optional): Months. If None, all months are estimated.
            xs (slice): X coordinates. If None, all x coordinates in source are estimated.
            ys (slice): Y coordinates. If None, all y coordinates in source are estimated.
            use_real_data (bool, optional): If available, use real data for estimation. Defaults to False.

        Returns:
            xr.Dataset: Dataset with wind speed.
        """

    @abc.abstractmethod
    def _estimate_cutout(self, *args, **kwargs) -> xr.Dataset:
        """Estimate the wind speed from a cutout.

        Args:
            height (int): Height of the wind speed, need to be greater than 0.
            years (slice): Years.
            months (slice, optional): Months. If None, all months are estimated.
            xs (slice): X coordinates. If None, all x coordinates in source are estimated.
            ys (slice): Y coordinates. If None, all y coordinates in source are estimated.
            use_real_data (bool, optional): If available, use real data for estimation. Defaults to False.

        Returns:
            xr.Dataset: Dataset with wind speed.
        """
