## Copyright 2020 Michael Davidson (UCSD), William Honaker.

## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License as
## published by the Free Software Foundation; either version 3 of the
## License, or (at your option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.

## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""
GEODATA

Geospatial Data Collection and "Pre-Analysis" Tools
"""

import os
import glob
import pandas as pd
import numpy as np
import xarray as xr
import shutil
import requests
from requests.exceptions import HTTPError
from six.moves import range
from contextlib import contextmanager
from tempfile import mkstemp
from calendar import monthrange

import logging
logger = logging.getLogger(__name__)

from ..config import modis_dir

datadir = modis_dir
spinup_var = False


def api_modis(
    toDownload,
    downloadedFiles
    ):
    if len(toDownload) == 0:
        logger.info("All MERRA2 files for this dataset have been downloaded.")
    else:
        count = 0
        print("get data goes here")


weather_data_config = {
	'modis_land_cover': dict(
		api_func=api_modis,
		file_granularity="yearly",
        band="LC_Type1",
		fn = os.path.join(modis_dir, '{year}/MODIS_006_MCD12Q1_{year}_01_01.nc')
	)
}