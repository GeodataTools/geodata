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
import hashlib
import json
import shutil
from pathlib import Path
from typing import Optional, Union

import xarray as xr

import geodata

from ..config import model_dir
from ..cutout import Cutout
from ..dataset import Dataset
from ..logging import logger
from ..utils import NpEncoder


class BaseModel(abc.ABC):
    """Base class for geospatial modeling.

    Args:
        name (str): The name of the model.
        source (Union[geodata.Dataset, geodata.Cutout, xr.Dataset]): The source of the model. Can be a Dataset, Cutout or xarray.Dataset.
        interpolate (bool, optional): Interpolate the source to the same grid as the target. Defaults to False.
        **kwargs: Additional keyword arguments to pass to the model.
    """

    SUPPORTED_WEATHER_DATA_CONFIGS: tuple[str]
    metadata_keys: set[str] = {
        "name",
        "module",
        "from_dataset",
        "years",
        "months",
        "files_orig",
        "files_prepared",
        "weather_data_config",
    }

    def __init__(self, source: Union[geodata.Dataset, geodata.Cutout], **kwargs):
        if source.config not in self.SUPPORTED_WEATHER_DATA_CONFIGS:
            raise ValueError(
                f"Weather data config {source.config} is not supported by this model."
            )

        if not source.prepared:
            raise ValueError(
                "The source Dataset/Cutout for this model is not prepared."
            )

        self.source = source
        self._extra_kwargs = kwargs
        self._corrupt_metadata = False
        self._prepared = False

        if isinstance(source, geodata.Dataset):
            self._ref_path = model_dir.parent / self.source.module
            self._path = model_dir / self.type / self.source.module
        else:
            self._ref_path = model_dir.parent / "cutouts"
            self._path = model_dir / self.type / self.source.name

        if (meta_path := self._path / "meta.json").exists():
            try:
                with open(meta_path, encoding="utf_8") as f:
                    self.metadata = json.load(f)
            except json.JSONDecodeError:
                logger.warning(
                    "Metadata file %s is corrupted. Model will be re-prepared.",
                    meta_path,
                )
                self._corrupt_metadata = True
                if isinstance(source, geodata.Dataset):
                    self.metadata = self.extract_dataset_metadata(source)
                else:
                    self.metadata = self.extract_cutout_metadata(source)
        else:
            # NOTE: We only store metadata in transient fashion until prepartion is done
            if isinstance(source, geodata.Dataset):
                self.metadata = self.extract_dataset_metadata(source)
            else:
                self.metadata = self.extract_cutout_metadata(source)

    def __repr__(self):
        return f"Model(source={self.source}, type={self.type})"

    @property
    def name(self) -> str:
        """Name of the model."""
        return f"{self.metadata['name']}_{self.type}"

    @property
    def module(self) -> str:
        """Module of the model."""
        return self.metadata["module"]

    @property
    @abc.abstractmethod
    def type(self) -> str:
        """Type of the model."""

    def estimate(
        self,
        height: int,
        years: slice,
        months: Optional[slice] = None,
        xs: Optional[slice] = None,
        ys: Optional[slice] = None,
        use_real_data: bool = False,
    ) -> xr.DataArray:
        """Estimate the wind speed at given coordinates.

        Args:
            height (int): Height of the wind speed, need to be greater than 0.
            years (slice): Years.
            months (slice, optional): Months. If None, all months are estimated.
            xs (slice): X coordinates. If None, all x coordinates in source are estimated.
            ys (slice): Y coordinates. If None, all y coordinates in source are estimated.
            use_real_data (bool, optional): If available, use real data for estimation. Defaults to False.

        Returns:
            xr.DataArray: Dataset with wind speed.
        """

        if self.from_dataset:
            return self._estimate_dataset(
                height=height,
                years=years,
                months=months,
                xs=xs,
                ys=ys,
                use_real_data=use_real_data,
            )

        return self._estimate_cutout(
            height=height,
            years=years,
            months=months,
            xs=xs,
            ys=ys,
            use_real_data=use_real_data,
        )

    def extract_dataset_metadata(self, dataset: Dataset) -> dict:
        if not dataset.prepared:
            raise ValueError("The source dataset for this model is not prepared.")

        logger.info("Using dataset %s", dataset.module)

        metadata = {}

        metadata["name"] = metadata["module"] = dataset.module
        metadata["from_dataset"] = True

        if isinstance(dataset.years, slice):
            metadata["years"] = dataset.years.start, dataset.years.stop
        if isinstance(dataset.months, slice):
            metadata["months"] = dataset.months.start, dataset.months.stop

        metadata["weather_data_config"] = dataset.config

        # NOTE: file paths for estimation parameters will be added later in the prepare step
        metadata["files_prepared"] = {}
        metadata["files_orig"] = {}
        for c, fp in dataset.downloadedFiles:
            if c == metadata["weather_data_config"]:
                with open(fp, "rb") as f:
                    metadata["files_orig"][
                        str(Path(fp).relative_to(self._ref_path))
                    ] = hashlib.sha256(f.read()).hexdigest()

        return metadata

    def extract_cutout_metadata(self, cutout: Cutout) -> dict:
        logger.info("Using cutout %s", cutout.name)

        metadata = {}

        metadata["name"] = cutout.name
        metadata["module"] = cutout.meta.attrs["module"]
        metadata["from_dataset"] = False
        metadata["weather_data_config"] = cutout.config

        if isinstance(cutout.years, slice):
            metadata["years"] = cutout.years.start, cutout.years.stop
        if isinstance(cutout.months, slice):
            metadata["months"] = cutout.months.start, cutout.months.stop

        metadata["files_prepared"] = {}
        metadata["files_orig"] = {}
        for yearmonth in cutout.coords["year-month"].to_index():
            fp = cutout.datasetfn(yearmonth)
            with open(fp, "rb") as f:
                metadata["files_orig"][str(Path(fp).relative_to(self._ref_path))] = (
                    hashlib.sha256(f.read()).hexdigest()
                )

        return metadata

    @property
    def prepared(self) -> bool:
        """Check if the model is prepared.

        Returns:
            bool: True if prepared.
        """

        if not self._prepared:
            self._prepared = self._check_prepared()
        return self._prepared

    def _check_prepared(self) -> bool:
        assert self._path is not None, "The model saving path has not been set yet."

        nc4_path: Path = self._path / "nc4"
        meta_path: Path = self._path / "meta.json"

        if not nc4_path.exists() or not meta_path.exists() or self._corrupt_metadata:
            return False
        with open(meta_path, encoding="utf-8") as f:
            metadata_loaded = json.load(f)
            if set(metadata_loaded.keys()) != self.metadata_keys:
                return False

        if len(self.metadata["files_orig"]) != len(self.metadata["files_prepared"]):
            return False

        nc4_rel_path = nc4_path.relative_to(self._path)
        for fp in self.metadata["files_orig"]:
            fp_prepared = str(nc4_rel_path / Path(fp).with_suffix(".params.nc4"))

            if (
                fp not in self.metadata["files_orig"]
                or fp_prepared not in self.metadata["files_prepared"]
            ):
                return False

            with open(self._ref_path / fp, "rb") as f:
                if (
                    self.metadata["files_orig"][fp]
                    != hashlib.sha256(f.read()).hexdigest()
                ):
                    logger.warning(
                        "File %s in source dataset has been modified since model creation. Model is not prepared!",
                        fp,
                    )
                    return False

            with open(self._path / fp_prepared, "rb") as f:
                if (
                    self.metadata["files_prepared"][fp_prepared]
                    != hashlib.sha256(f.read()).hexdigest()
                ):
                    logger.warning(
                        "Parameter file %s in model has been modified since model creation. Model is not prepared!",
                        fp_prepared,
                    )
                    return False

        return True

    def prepare(self, force: bool = False):
        """Prepare the model.

        Args:
            force (bool, optional): Force re-prepare the model. Defaults to False."""

        self._prepared = False  # NOTE: force re-checking preparedness here
        if self.prepared and not force:
            logger.info("The model is already prepared.")
            return

        logger.info("Model not present in model directory, creating.")
        shutil.rmtree(self._path, ignore_errors=True)
        (self._path / "nc4").mkdir(exist_ok=True, parents=True)

        self.metadata["files_prepared"] = {}
        for fp in (
            self._prepare_dataset() if self.from_dataset else self._prepare_cutout()
        ):
            with open(self._path / fp, "rb") as f:
                self.metadata["files_prepared"][fp] = hashlib.sha256(
                    f.read()
                ).hexdigest()

        with open(self._path / "meta.json", "w", encoding="utf-8") as f:
            json.dump(self.metadata, f, indent=4, cls=NpEncoder)

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
        if "from_dataset" not in self.metadata:
            return False
        return self.metadata["from_dataset"]

    @abc.abstractmethod
    def _estimate_dataset(
        self,
        height: int,
        years: slice,
        months: Optional[slice] = None,
        xs: Optional[slice] = None,
        ys: Optional[slice] = None,
        use_real_data: Optional[bool] = False,
    ) -> xr.DataArray:
        """Estimate the wind speed from a dataset.

        Args:
            height (int): Height of the wind speed, need to be greater than 0.
            years (slice): Years.
            months (slice, optional): Months. If None, all months are estimated.
            xs (slice): X coordinates. If None, all x coordinates in source are estimated.
            ys (slice): Y coordinates. If None, all y coordinates in source are estimated.
            use_real_data (bool, optional): If available, use real data for estimation. Defaults to False.

        Returns:
            xr.DataArray: Dataset with wind speed.
        """

    @abc.abstractmethod
    def _estimate_cutout(
        self,
        height: int,
        years: slice,
        months: Optional[slice] = None,
        xs: Optional[slice] = None,
        ys: Optional[slice] = None,
        use_real_data: Optional[bool] = False,
    ) -> xr.DataArray:
        """Estimate the wind speed from a cutout.

        Args:
            height (int): Height of the wind speed, need to be greater than 0.
            years (slice): Years.
            months (slice, optional): Months. If None, all months are estimated.
            xs (slice): X coordinates. If None, all x coordinates in source are estimated.
            ys (slice): Y coordinates. If None, all y coordinates in source are estimated.
            use_real_data (bool, optional): If available, use real data for estimation. Defaults to False.

        Returns:
            xr.DataArray: Dataset with wind speed.
        """

    @property
    def files(self):
        if "files_prepared" not in self.metadata or "files_orig" not in self.metadata:
            return []

        files_prepared = [self._path / p for p in self.metadata["files_prepared"]]
        files_orig = [self._ref_path / p for p in self.metadata["files_orig"]]

        return files_prepared + files_orig
