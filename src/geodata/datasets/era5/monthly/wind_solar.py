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
import os
import pprint
import tempfile
import zipfile
from pathlib import Path

import xarray as xr

from ..hourly.wind_solar import ERA5WindSolarHourlyDataset

logger = logging.getLogger(__name__)


class ERA5WindSolarMonthlyDataset(ERA5WindSolarHourlyDataset):
    """ERA5WindSolarMonthlyDataset is a class that handles the downloading,
    preprocessing, and storing of the ERA5 dataset for wind and solar
    information. This dataset is stored in monthly intervals.

    The ERA5 dataset is a reanalysis dataset that provides a comprehensive
    record of the Earth's climate. It is produced by the European Centre for
    Medium-Range Weather Forecasts (ECMWF) and is available from 1980 to
    present.

    Note:
        - The specific variables that are downloaded are:
            - 100m_u_component_of_wind
            - 100m_v_component_of_wind
            - 2m_temperature
            - runoff
            - soil_temperature_level_4
            - surface_net_solar_radiation
            - surface_pressure
            - surface_solar_radiation_downwards
            - toa_incident_solar_radiation
            - total_sky_direct_solar_radiation_at_surface
            - forecast_surface_roughness
            - geopotential
    """

    weather_config = "wind_solar_monthly"

    def _download_file(self, file: dict):
        year: int = file["year"]
        month: int = file["month"]
        save_path: Path = file["save_path"]

        full_request = {
            "product_type": self.product_type,
            "format": "netcdf",
            "variable": list(self.variables.keys()),
            "year": year,
            "month": month,
            "day": [f"{d:02d}" for d in range(1, 32)],
            "time": "00:00",
        }

        if self.bounds is not None:
            full_request["area"] = self.bounds[::-1]

        logger.debug("Full request for download: %s", pprint.pformat(full_request))

        full_result = self.client.retrieve(self.product, full_request)
        if full_result.content_type == "application/zip":
            logger.info(
                "Multiple files found with request. Additional unzipping/preprocessing needed."
            )

            with tempfile.TemporaryDirectory() as tempdir:
                full_result.download(os.path.join(tempdir, "download.zip"))
                with zipfile.ZipFile(
                    os.path.join(tempdir, "download.zip"), "r"
                ) as zip_ref:
                    zip_ref.extractall(tempdir)

                with xr.open_mfdataset(
                    [
                        os.path.join(tempdir, f)
                        for f in os.listdir(tempdir)
                        if f.endswith(".nc")
                    ]
                ) as ds:
                    ds.to_netcdf(save_path)

                logger.info("Preprocessing complete with zipfile")
                logger.info("Successfully downloaded to %s", save_path)
