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


import hashlib
import logging
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path

import xarray as xr
from tqdm.auto import tqdm

from geodata.utils import check_hash

from ._base import BaseModelResult

logger = logging.getLogger(__name__)


@dataclass
class DailyModelResult(BaseModelResult):
    """Class for daily model results."""

    def _check_prepared(self):
        assert self.path is not None, "The model saving path has not been set yet."

        if not (self.path / "meta.json").exists():
            logger.warning(
                "Model %s-%s does not have metadata. Please prepare the model first.",
                self.year,
                self.month,
            )
            return False

        if self.model.quick_check:
            # If the quick check is enabled, we only need to check the file hashes
            # and not the actual data.
            logger.debug(
                "Quick check is enabled. Only checking file hashes for model %s-%s.",
                self.year,
                self.month,
            )

            for file in self.files:
                if not file.exists():
                    logger.warning(
                        "File %s in model does not exist. Model is not prepared!",
                        str(file),
                    )
                    return False

            return True
        from .._base import MAX_WORKERS

        with ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
            files = [(f, self._hashes.get(f.name)) for f in self.files]
            results = list(
                tqdm(
                    executor.map(lambda t: check_hash(*t), files),
                    total=len(files),
                    unit="file",
                    dynamic_ncols=True,
                    desc=f"Model Files Integrity Check {self.year:04d}-{self.month:02d}",
                )
            )

        for file, (is_valid, hash_value) in zip(files, results):
            if not is_valid:
                logger.warning(
                    "File %s in model has been modified since model creation. Model is not prepared!",
                    str(file[0]),
                )
                return False

        return True

    def register(self, dataset: xr.Dataset):
        """Register the model result with the dataset. This file should be a file
        covering the entire model period.

        Args:
            dataset (xr.Dataset): The dataset to register with.
        """

        time = dataset.get("valid_time")

        start = time.min()
        end = time.max()

        if start.dt.month != self.month or end.dt.month != self.month:
            raise ValueError(
                f"Dataset start month {start.dt.month} does not match model month {self.month}."
            )

        if start.dt.year != self.year or end.dt.year != self.year:
            raise ValueError(
                f"Dataset start year {start.dt.year} does not match model year {self.year}."
            )

        days, datasets = zip(*dataset.groupby("valid_time.day"))
        paths = [self.path / f"{day:02d}.nc" for day in days]

        logger.debug("Saving model results to %s", self.path)

        from .._base import XR_ENGINE

        xr.save_mfdataset(datasets, paths, engine=XR_ENGINE)

        # Write the hash file for integrity checking
        with ThreadPoolExecutor() as executor:

            def compute_hash(path: Path):
                sha256 = hashlib.sha256()
                with path.open("rb") as f:
                    for chunk in iter(lambda: f.read(8192), b""):
                        sha256.update(chunk)
                return path.name, sha256.hexdigest()

            results = list(
                tqdm(
                    executor.map(compute_hash, paths),
                    total=len(paths),
                    unit="file",
                    dynamic_ncols=True,
                    desc="Computing File Hashes",
                )
            )

        for k, v in results:
            self._hashes[k] = v
