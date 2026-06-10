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
from dataclasses import dataclass

import xarray as xr

from geodata.utils import check_hash

from ._base import BaseModelResult

logger = logging.getLogger(__name__)


@dataclass
class MonthlyModelResult(BaseModelResult):
    """Class for monthly model results."""

    def _check_prepared(self):
        assert self.path is not None, "The model saving path has not been set yet."

        if not (self.path / "meta.json").exists():
            logger.warning(
                "Model %s-%s does not have metadata. Please prepare the model first.",
                self.year,
                self.month,
            )
            return False

        return check_hash(self.path / f"{self.month:02d}.params.nc")[0]

    def register(self, dataset: xr.Dataset):
        from .._base import _get_xr_engine
        
        engine = _get_xr_engine()
        logger.info(f"register: Saving monthly file with engine={engine}")
        dataset.to_netcdf(self.path / f"{self.month:02d}.params.nc", engine=engine)
        with open(self.path / f"{self.month:02d}.params.nc", "rb") as f:
            self._hashes[f"{self.month:02d}.params.nc"] = hashlib.sha256(
                f.read()
            ).hexdigest()
