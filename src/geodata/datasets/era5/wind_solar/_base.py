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

import xarray as xr
import numpy as np

from geodata.types import PathLike
from .._base import ERA5BaseDataset, _subset_x_y_era5

logger = logging.getLogger(__name__)

def _add_height(ds):
    """Convert geopotential 'z' to geopotential height following [1]

    References
    ----------
    [1] ERA5: surface elevation and orography, retrieved: 10.02.2019
    https://confluence.ecmwf.int/display/CKB/ERA5%3A+surface+elevation+and+orography

    """
    g0 = 9.80665
    z = ds["z"]
    if "time" in z.coords:
        z = z.isel(time=0, drop=True)
    ds["height"] = z / g0
    ds = ds.drop("z")
    return ds

class ERA5WindSolarBaseDataset(ERA5BaseDataset):
    """Base class for ERA5 wind and solar datasets.
    
    This class provides the prepare_func implementation specific to wind_solar datasets,
    which use single-level data from the reanalysis-era5-single-levels product.
    """

    @classmethod
    def transform_wind_solar_dataset(cls, ds: xr.Dataset) -> xr.Dataset:
        """Transform raw ERA5 wind_solar dataset to standardized variable names and units.
        
        This method applies the transformations needed to convert raw ERA5 variables
        to the standardized format used by models (e.g., pvlib). It can be called
        directly on an already-opened dataset.
        
        Args:
            ds: Raw ERA5 dataset with original variable names (fdir, tisr, t2m, u100, v100, etc.)
        
        Returns:
            Transformed dataset with standardized variable names (influx_direct, influx_diffuse,
            temperature, wnd100m, etc.)
        """
        # Add height from geopotential if not already present
        if "height" not in ds.data_vars and "z" in ds.data_vars:
            ds = _add_height(ds)
        
        # Rename radiation variables
        ds = ds.rename({"fdir": "influx_direct", "tisr": "influx_toa"})
        
        # Calculate albedo and influx_diffuse
        with np.errstate(divide="ignore", invalid="ignore"):
            ds["albedo"] = (
                ((ds["ssrd"] - ds["ssr"]) / ds["ssrd"])
                .fillna(0.0)
                .assign_attrs(units="(0 - 1)", long_name="Albedo")
            )
        influx_diffuse = ds["ssrd"] - ds["influx_direct"]
        influx_diffuse.attrs.update({
            "units": "J m**-2",
            "long_name": "Surface diffuse solar radiation downwards"
        })
        ds["influx_diffuse"] = influx_diffuse
        ds = ds.drop(["ssrd", "ssr"])

        # Convert from energy to power J m**-2 -> W m**-2 and clip negative fluxes
        for a in ("influx_direct", "influx_diffuse", "influx_toa"):
            ds[a] = ds[a].clip(min=0.0) / (60.0 * 60.0)
            ds[a].attrs["units"] = "W m**-2"

        # Calculate wind speed from u and v components
        wnd100m = np.sqrt(ds["u100"] ** 2 + ds["v100"] ** 2)
        if isinstance(wnd100m, xr.DataArray):
            wnd100m.attrs.update({
                "units": ds["u100"].attrs.get("units", ""),
                "long_name": "100 metre wind speed"
            })
        else:
            # If it's a numpy array, convert to DataArray with attrs
            wnd100m = xr.DataArray(
                wnd100m,
                coords=ds["u100"].coords,
                dims=ds["u100"].dims,
                attrs={
                    "units": ds["u100"].attrs.get("units", ""),
                    "long_name": "100 metre wind speed"
                }
            )
        ds["wnd100m"] = wnd100m
        ds = ds.drop(["u100", "v100"])

        # Rename other variables
        ds = ds.rename(
            {
                "ro": "runoff",
                "t2m": "temperature",
                "sp": "pressure",
                "stl4": "soil temperature",
                "fsr": "roughness",
                "d2m": "dewpoint_temperature"
            }
        )

        # New ERA5 format for hourly datasets
        # See https://forum.ecmwf.int/t/new-time-format-in-era5-netcdf-files/3796
        # TODO: We can remove this if we refactor geodata's convert module in the future
        if "valid_time" in ds.coords:
            ds = ds.rename({"valid_time": "time"})
        
        return ds

    @classmethod
    def prepare_func(
        cls,
        fn: PathLike,
        year: int,
        month: int,
        xs: slice,
        ys: slice,
        **kwargs,
    ):
        """Prepare the dataset for a given year and month.
        
        This implementation is specific to wind_solar datasets which:
        - Use single-level data (no model levels)
        - Download from reanalysis-era5-single-levels product
        - Are stored as monthly files
        """
        if isinstance(fn, str) and not os.path.exists(fn):
            return
        if isinstance(fn, list) and not all(os.path.isfile(f) for f in fn):
            return

        with xr.open_dataset(fn) as ds:
            logger.info("Opening %s", fn)
            ds = _add_height(ds)
            ds = _subset_x_y_era5(ds, xs, ys)
            
            # Use the shared transformation method
            ds = cls.transform_wind_solar_dataset(ds)

            yield (year, month), ds

