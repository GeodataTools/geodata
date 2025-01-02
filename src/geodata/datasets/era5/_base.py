# Copyright 2024-2025 Michael Davidson (UCSD), Xiqiang Liu (UCSD)

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

import logging

import cdsapi

from .._base import BaseDataset


class ERA5BaseDataset(BaseDataset):
    """ERA5BaseDataset is a class that encaps a dataset from the ERA5 reanalysis
    dataset. It provides a streamlined workflow for downloading, preprocessing,
    and storing of these datasets.
    """

    module = "era5"

    def _extra_setup(self, **kwargs):
        self.logger = logging.getLogger(__name__.replace("._base", ".client"))
        self.client = cdsapi.Client(
            info_callback=self.logger.info,
            error_callback=self.logger.error,
            debug_callback=self.logger.debug,
            warning_callback=self.logger.warning,
        )
