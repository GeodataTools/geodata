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
from pathlib import Path

import xarray as xr

from ..._base import AtomicDataset
from ._base import ERA5Wind3DBaseDataset

logger = logging.getLogger(__name__)


class ERA5Wind3DHourlyDataset(ERA5Wind3DBaseDataset):
    """ERA5Wind3DHourlyDataset is a class that handles the downloading,
    preprocessing, and storing of the ERA5 dataset for wind information.
    This dataset is stored in hourly intervals.
    """

    weather_config = "wind_3d_hourly"
    frequency = "daily"
    L137_LEVELS = range(131, 138)

    # Information that are needed for ERA5's API request
    variables = {"u": "131", "v": "132"}
    product = "reanalysis-era5-complete"

    def _download_file(self, file: AtomicDataset):
        """Sample request for reference:

        c.retrieve('reanalysis-era5-complete', {
        'date'    : '20130101/to/20130131',
        'levelist': '1/10/100/137',
        'levtype' : 'ml',
        'param'   : '130,  # Full information at https://apps.ecmwf.int/codes/grib/param-db/
        'stream'  : 'oper',  # Denotes ERA5. Ensemble members are selected by 'enda'
        'time'    : '00/to/23/by/6',
        'type'    : 'an',
        'grid'    : '1.0/1.0',
        'format'  : 'netcdf',
        }, 'save_path.nc')  # Output file. Adapt as you wish.
        """

        year: int = file.year
        month: int = file.month
        day: int = file.day
        save_path: Path = file.path

        full_request = {
            "date": f"{year}{month:02d}{day:02d}",
            "levelist": "/".join([str(level) for level in self.L137_LEVELS]),
            "levtype": "ml",
            "param": "/".join(self.variables.values()),
            "stream": "oper",
            "time": "00/to/23/by/1",
            "type": "an",
            "grid": "0.25/0.25",  # NOTE: We want the highest resolution possible
            "format": "netcdf",
        }

        logger.debug("Full request for download: %s", pprint.pformat(full_request))

        full_result = self.client.retrieve(self.product, full_request)

        if not self.bounds:
            _count = 0
            for _ in range(3):
                try:
                    _count += 1
                    full_result.download(save_path)
                    # Ensure file is fully written to disk before proceeding
                    with open(save_path, "rb") as f:
                        os.fsync(f.fileno())
                    logger.info("File downloaded: %s", save_path)
                    return
                except Exception as e:
                    logger.error("Download failed: %s", e)
                    if _count == 3:
                        raise

        # NOTE: Raw MARS request doesn't support bounding box, so we need to
        # subset the data after downloading
        with tempfile.NamedTemporaryFile(suffix=".nc") as tmpfile:
            _count = 0
            for _ in range(3):
                try:
                    _count += 1
                    full_result.download(tmpfile.name)
                    # Ensure file is fully written to disk before proceeding
                    os.fsync(tmpfile.fileno())
                    logger.info("File downloaded: %s", save_path)
                    break
                except Exception as e:
                    logger.error("Download failed: %s", e)
                    if _count == 3:
                        raise
            with xr.open_dataset(tmpfile.name, chunks="auto", engine="h5netcdf") as ds:
                ds = ds.sel(
                    longitude=slice(*sorted([self.bounds[0], self.bounds[2]])),
                    latitude=slice(
                        *sorted([self.bounds[1], self.bounds[3]], reverse=True)
                    ),
                )
                ds.to_netcdf(save_path, engine="h5netcdf")
                # Ensure file is fully written to disk before proceeding
                with open(save_path, "rb") as f:
                    os.fsync(f.fileno())

        logger.info("File downloaded: %s", save_path)

