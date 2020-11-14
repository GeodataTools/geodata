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
import ee
import oauth2client
import numpy as np
from osgeo import gdal
from osgeo import osr
import time
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
    # ee.Authenticate()
    ee.Initialize()
    if len(toDownload) == 0:
        logger.info("All MERRA2 files for this dataset have been downloaded.")
    else:
        count = 0
        
        for f in toDownload:
            nc_filename = f[1]
            tif_filename = f[2]
            img = ee.Image(f[3])
            area = ee.Geometry.Rectangle(f[4])
            band = f[5]

            img = img.select(band)

            latlon = ee.Image.pixelLonLat().addBands(img)
            latlon_r = latlon.reduceRegion(
                reducer=ee.Reducer.toList(),
                geometry=area,
                maxPixels=1e9,
                scale=1000)

            data = np.array((ee.Array(latlon_r.get(band)).getInfo()))
            lats = np.array((ee.Array(latlon_r.get("latitude")).getInfo()))
            lons = np.array((ee.Array(latlon_r.get("longitude")).getInfo()))

            uniqueLats = np.unique(lats)
            uniqueLons = np.unique(lons)
            ncols = len(uniqueLons)    
            nrows = len(uniqueLats)

            ys = uniqueLats[1] - uniqueLats[0] 
            xs = uniqueLons[1] - uniqueLons[0]

            arr = np.zeros([nrows, ncols], np.float32) 

            counter =0
            for y in range(0,len(arr),1):
                for x in range(0,len(arr[0]),1):
                    if lats[counter] == uniqueLats[y] and lons[counter] == uniqueLons[x] and counter < len(lats)-1:
                        counter+=1
                        arr[len(uniqueLats)-1-y,x] = data[counter] 

            transform = (np.min(uniqueLons),xs,0,np.max(uniqueLats),0,-ys)

            target = osr.SpatialReference()
            target.ImportFromEPSG(4326)

            driver = gdal.GetDriverByName('GTiff')
            timestring = time.strftime("%Y%m%d_%H%M%S")

            outputDataset = driver.Create(tif_filename, ncols, nrows, 1, gdal.GDT_Float32)
            outputDataset.SetMetadata( {
                'filename': nc_filename,
                'image': f[3],
                'band': band,
                'area': f[4]
                })
            outputDataset.SetGeoTransform(transform)
            outputDataset.SetProjection(target.ExportToWkt())
            outputDataset.GetRasterBand(1).WriteArray(arr)
            outputDataset.GetRasterBand(1).SetNoDataValue(-9999)
            outputDataset = None

            gdal_ds = gdal.Translate(nc_filename, tif_filename, format='NetCDF')
            
            downloadedFiles.append((f[0], nc_filename))
            count += 1

            logger.info("Successfully downloaded data for %s", nc_filename)
        #


weather_data_config = {
	'modis_land_cover': dict(
		api_func=api_modis,
		file_granularity="yearly",
        band="LC_Type1",
        image='MODIS/006/MCD12Q1/{year}_01_01',
		nc_fn = os.path.join(modis_dir, '{year}/MODIS_006_MCD12Q1_{year}_01_01.nc'),
        tif_fn = os.path.join(modis_dir, '{year}/MODIS_006_MCD12Q1_{year}_01_01.nc')
	)
}