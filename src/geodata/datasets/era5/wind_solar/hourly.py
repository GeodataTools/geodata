# Copyright 2024-2025 Michael Davidson (UCSD), Xiqiang Liu (UCSD), Keyu Long (UCSD)

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

from ..._base import AtomicDataset
from ._base import ERA5WindSolarBaseDataset

logger = logging.getLogger(__name__)


class ERA5WindSolarHourlyDataset(ERA5WindSolarBaseDataset):
    """ERA5WindSolarHourlyDataset is a class that handles the downloading,
    preprocessing, and storing of the ERA5 dataset for wind and solar
    information. This dataset is stored in hourly intervals.

    The ERA5 dataset is a reanalysis dataset that provides a comprehensive
    record of the Earth's climate. It is produced by the European Centre for
    Medium-Range Weather Forecasts (ECMWF) and is available from 1980 to
    present.

    Note:
        - The specific variables that are downloaded are:
            - 100m_u_component_of_wind
            - 100m_v_component_of_wind
            - 2m_temperature
            - 2m_dewpoint_temperature
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

    weather_config = "wind_solar_hourly"

    # Information that is needed for ERA5's API request
    variables = {
        "100m_u_component_of_wind": "u100",
        "100m_v_component_of_wind": "v100",
        "2m_temperature": "t2m",
        "2m_dewpoint_temperature": "d2m",
        "runoff": "ro",
        "soil_temperature_level_4": "stl4",
        "surface_net_solar_radiation": "ssr",
        "surface_pressure": "sp",
        "surface_solar_radiation_downwards": "ssrd",
        "toa_incident_solar_radiation": "tisr",
        "total_sky_direct_solar_radiation_at_surface": "fdir",
        "forecast_surface_roughness": "fsr",
        "geopotential": "z",
    }
    product = "reanalysis-era5-single-levels"
    product_type = "reanalysis"

    def _download_file(self, file: AtomicDataset):
        year: int = file.year
        month: int = file.month
        save_path: Path = file.path

        # Limit to first 3 days when testing=True
        max_day = 4 if self.testing else 32
        
        full_request = {
            "product_type": self.product_type,
            "format": "netcdf",
            "variable": list(self.variables.keys()),
            "year": year,
            "month": month,
            "day": [f"{d:02d}" for d in range(1, max_day)],
            "time": [f"{t:02d}:00" for t in range(0, 24)],
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
                _count = 0
                for _ in range(3):
                    try:
                        _count += 1
                        full_result.download(os.path.join(tempdir, "download.zip"))
                        break
                    except Exception as e:
                        logger.error("Error downloading file: %s", e)
                        if _count == 3:
                            raise

                with zipfile.ZipFile(
                    os.path.join(tempdir, "download.zip"), "r"
                ) as zip_ref:
                    zip_ref.extractall(tempdir)

                with xr.open_mfdataset(
                    [
                        os.path.join(tempdir, f)
                        for f in os.listdir(tempdir)
                        if f.endswith(".nc")
                    ], engine="h5netcdf"
                ) as ds:
                    ds.to_netcdf(save_path, engine="h5netcdf")

                logger.info("Preprocessing complete with zipfile")
                logger.info("Successfully downloaded to %s", save_path)

