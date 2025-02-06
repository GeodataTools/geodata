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

import pandas as pd
import xarray as xr
from herbie import FastHerbie, Herbie

from ...logging import redirect_stdout_to_logger
from .._base import AtomicDataset
from ._base import HRRRBaseDataset

logger = logging.getLogger(__name__)


class HRRRHourlyDataset(HRRRBaseDataset):
    """HRRRHourlyDataset is a class that encaps a dataset from the HRRR
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

    weather_config = "wind_solar"
    product = "sfc"

    def _download_file(self, file: AtomicDataset):
        year, month, day = file.year, file.month, file.day

        date_range = pd.date_range(
            f"{year}-{month}-{day}",
            f"{year}-{month}-{day+1}",
            freq="h",
            inclusive="left",
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
            fh.download(":[UV]GRD:[1,8]0 m")
            fh.download(":TMP:2 m")

            uv_10 = []
            uv_80 = []
            tmp_2 = []
            wrfs = []

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
                        self._preprocess_individual_herbie(h, "[UV]GRD:10 m", hour)
                    )
                    uv_80.append(
                        self._preprocess_individual_herbie(h, "[UV]GRD:80 m", hour)
                    )
                    tmp_2.append(
                        self._preprocess_individual_herbie(h, ":TMP:2 m", hour)
                    )
                    wrfs.append(
                        self._preprocess_individual_herbie(h, ":..WRF:surface", hour)
                    )

                except ValueError:
                    logger.warning(f"No data found for {hour}, skipping.")

            # First concat to daily files

            uv_10: xr.Dataset = xr.open_mfdataset(
                uv_10, concat_dim="time", chunks="auto", combine="nested"
            )
            uv_80: xr.Dataset = xr.open_mfdataset(
                uv_80, concat_dim="time", chunks="auto", combine="nested"
            ).rename({"u": "u80", "v": "v80"})
            tmp_2: xr.Dataset = xr.open_mfdataset(
                tmp_2, concat_dim="time", chunks="auto", combine="nested"
            )
            wrfs: xr.Dataset = xr.open_mfdataset(
                wrfs, concat_dim="time", chunks="auto", combine="nested"
            )

            ds = xr.merge([uv_10, uv_80, tmp_2, wrfs], compat="override")

            try:
                del ds["heightAboveGround"]
                del ds.attrs["search"]
                del ds.attrs["local_grib"]
                del ds.attrs["remote_grib"]
                del ds.coords["gribfile_projection"]
            except KeyError:
                pass

            ds.to_netcdf(file.path)

        # NOTE: Flush temporary FastHerbie save directory to save space, since we no
        # longer need the raw downloaded files
        self._herbie_save_dir.cleanup()
