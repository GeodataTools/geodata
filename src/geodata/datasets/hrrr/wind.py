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

import logging
import multiprocessing as mp
import xarray as xr

import pandas as pd
from herbie import FastHerbie, Herbie

from ...logging import redirect_stdout_to_logger
from ._base import HRRRBaseDataset

logger = logging.getLogger(__name__)


class HRRRWindHourlyDataset(HRRRBaseDataset):
    """HRRRWindHourlyDataset is a class that encaps a dataset from the HRRR
    dataset. It provides a streamlined workflow for downloading, preprocessing,
    and storing of these datasets.

    The HRRR dataset is a high-resolution weather forecast model that provides
    hourly data for the United States. This dataset is useful for a variety of
    applications, including renewable energy forecasting, weather prediction,
    and climate research.

    This class provides a simple interface for downloading and processing the
    HRRR dataset. It allows users to specify the years and months of interest,
    as well as the variables they wish to download.
    """

    weather_config = "wind"
    product = "sfc"

    def _download_file(self, file: dict):
        year, month = file["year"], file["month"]

        date_range = pd.date_range(
            f"{year}-{month}-01", f"{year}-{month+1}-01", freq="h", inclusive="left"
        )

        fh = FastHerbie(
            date_range,
            model=self.module,
            product=self.product,
            max_threads=mp.cpu_count() * 2,
            save_dir=self._herbie_save_dir.name,
            priority=self._priority,
        )

        with redirect_stdout_to_logger(logger, logging.INFO):
            logger.info(f"Downloading HRRR wind data in bulk for {year}/{month}")
            fh.download(":[UV]GRD:[1,8]0 m")

            uv_10 = []
            uv_80 = []
            for hour in date_range:
                h = Herbie(
                    hour,
                    model=self.module,
                    product=self.product,
                    save_dir=self._herbie_save_dir.name,
                    priority=self._priority,
                )

                try:
                    uv_10.append(
                        h.xarray("[UV]GRD:10 m").rename({"u10": "u", "v10": "v"})
                    )
                    uv_80.append(h.xarray("[UV]GRD:80 m"))
                except ValueError:
                    logger.warning(f"No data found for {hour}, skipping.")

            uv_10 = xr.concat(uv_10, dim="time")
            uv_80 = xr.concat(uv_80, dim="time")

            ds: xr.Dataset = xr.concat([uv_10, uv_80], dim="heightAboveGround").rename(
                {"latitude": "y", "longitude": "x"}
            )
            ds["wind_speed"] = (ds["u"] ** 2 + ds["v"] ** 2) ** 0.5

            try:
                del ds.attrs["search"]
                del ds.attrs["local_grib"]
                del ds.attrs["remote_grib"]
            except KeyError:
                pass

            ds.to_netcdf(file["save_path"])

        # NOTE: Flush temporary FastHerbie save directory to save space, since we no
        # longer need the raw downloaded files
        self._herbie_save_dir.cleanup()
