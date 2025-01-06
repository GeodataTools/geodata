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

from pathlib import Path

import numpy as np
import xarray as xr
import requests

from .._base import BaseDataset


class MERRA2BaseDataset(BaseDataset):
    """MERRA2BaseDataset is a class that encaps a dataset from the MERRA2 reanalysis
    dataset. It provides a streamlined workflow for downloading, preprocessing,
    and storing of these datasets.
    """

    module = "merra2"
    url_template = ""

    def _download_file(self, file: dict):
        assert "url" in file, "URL is required to download the file"

        url: str = file["url"]
        path: Path = file["path"]

        # Download the file
        with requests.get(url, stream=True) as r:
            r.raise_for_status()

            with open(path, "wb") as f:
                for chunk in r.iter_content(chunk_size=8192):
                    f.write(chunk)

    def spinup_year(year: int, month: int):
        """Returns the spinup period for the given year and month.
        See https://gmao.gsfc.nasa.gov/pubs/docs/Bosilovich785.pdf for more
        information.

        Args:
            year (int): The year of the dataset
            month (int): The month of the dataset

        Returns:
            str: The spinup period
        """
        if year >= 1980 and year < 1992:
            spinup = "100"
        elif year >= 1992 and year < 2001:
            spinup = "200"
        elif year >= 2001 and year < 2011:
            spinup = "300"
        elif year >= 2011 and year < 2020:
            spinup = "400"
        elif year == 2020 and month == 9:
            spinup = "401"
        else:
            spinup = "400"

        return spinup

    def convert_and_subset_lons_lats_merra2(
        ds: xr.Dataset | xr.DataArray, xs: slice, ys: slice
    ):
        """Rename geographic dimensions to x,y. Subset x,y according to xs, ys.

        Args:
            ds (xr.Dataset | xr.DataArray): The dataset to subset
            xs (slice): The slice of longitudes to subset
            ys (slice): The slice of latitudes to subset

        Returns:
            xr.Dataset | xr.DataArray: The subsetted dataset
        """

        if not isinstance(xs, slice):
            first, second, last = np.asarray(xs)[[0, 1, -1]]
            xs = slice(first - 0.1 * (second - first), last + 0.1 * (second - first))
        if not isinstance(ys, slice):
            first, second, last = np.asarray(ys)[[0, 1, -1]]
            ys = slice(first - 0.1 * (second - first), last + 0.1 * (second - first))

        ds = ds.sel(lat=ys)

        # Longitudes should go from -180. to +180.
        if len(ds.coords["lon"].sel(lon=slice(xs.start + 360.0, xs.stop + 360.0))):
            ds = xr.concat(
                [ds.sel(lon=slice(xs.start + 360.0, xs.stop + 360.0)), ds.sel(lon=xs)],
                dim="lon",
            )
            ds = ds.assign_coords(
                lon=np.where(
                    ds.coords["lon"].values <= 180,
                    ds.coords["lon"].values,
                    ds.coords["lon"].values - 360.0,
                )
            )
        else:
            ds = ds.sel(lon=xs)

        ds = ds.rename({"lon": "x", "lat": "y"})
        ds = ds.assign_coords(lon=ds.coords["x"], lat=ds.coords["y"])

        return ds

    def _rename_and_clean_coords(
        ds: xr.Dataset | xr.DataArray, add_lon_lat: bool = True
    ):
        """Rename 'longitude' and 'latitude' columns to 'x' and 'y'

        Optionally `add_lon_lat` preserves latitude and longitude
        columns as 'lat' and 'lon'.

        Args:
            ds (xr.Dataset): The dataset to rename
            add_lon_lat (bool, optional): Whether to preserve latitude and longitude columns. Defaults to True.
        """

        ds = ds.rename({"lon": "x", "lat": "y"})
        if add_lon_lat:
            ds = ds.assign_coords(lon=ds.coords["x"], lat=ds.coords["y"])
        return ds

    def _hourly_catalog(self):
        if not self.url_template:
            raise NotImplementedError("url_template is not defined for this dataset")

        catalog = super()._hourly_catalog()

        for file in catalog:
            file["spinup"] = self.spinup_year(file["year"], file["month"])
            file["url"] = self.url_template.format(**file)

        return catalog
