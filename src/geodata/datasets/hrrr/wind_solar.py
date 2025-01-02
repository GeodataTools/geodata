# Copyright 2024 Michael Davidson (UCSD), Xiqiang Liu (UCSD)

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

from multiprocessing import cpu_count
from unittest.mock import patch


import pandas as pd
import xarray as xr
from herbie import Herbie
from tqdm.auto import tqdm

from ._base import HRRRBaseDataset

logger = logging.getLogger(__name__)


def fake_print(*args, **kwargs):
    pass


class HRRRWindSolarDataset(HRRRBaseDataset):
    """HRRRWindSolarDataset is a class that encaps a dataset from the HRRR
    dataset. It provides a streamlined workflow for downloading, preprocessing,
    and storing of these datasets.
    """

    weather_config = "wind_solar"
    variables = ":[UV]GRD:[1,8]0 m"

    def download(self):
        """Download the dataset from the HRRR dataset."""
        logger.info(f"Downloading {self.weather_config} dataset")
        logger.info(self._herbie_save_dir.name)

        with patch("builtins.print", fake_print):
            for file in tqdm(self.catalog, desc="Downloading Data", dynamic_ncols=True):
                hours = pd.date_range(
                    f"{file['year']}-{file['month']}-01",
                    f"{file['year']}-{file['month']}-31",
                    freq="1h",
                )

                dss = []
                for hour in hours:
                    h = Herbie(
                        hour,
                        fxx=0,
                        product="sfc",
                        model="hrrr",
                        priority=self._priority,
                        save_dir=self._herbie_save_dir.name,
                        max_threads=cpu_count() * 2,
                    )
                    logger.info(f"Downloading {hour}")
                    dss.append(
                        xr.concat(h.xarray(self.variables), dim="heightAboveGround")
                    )

                ds = xr.concat(dss, dim="time").to_netcdf(file["save_path"])
                ds.close()

        logger.info(f"Downloaded {self.weather_config} dataset")

    @property
    def downloaded(self):
        pass
