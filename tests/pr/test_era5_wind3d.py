# Copyright 2025 Keyu Long (UCSD)

# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License as
# published by the Free Software Foundation; either version 3 of the
# License, or (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import logging

import xarray as xr
from dask.distributed import Client

from geodata.datasets import load_dataset
from geodata.datasets._base import BaseDataset
from geodata.logging import logger
from geodata.model.wind import WindInterpolationModel

logger.setLevel(logging.DEBUG)


def _fixture_xy_slices(dataset: BaseDataset) -> tuple[slice, slice]:
    """Build ``xs``, ``ys`` slices on the fixture grid (``x``/``y`` or ERA5 ``longitude``/``latitude``)."""
    path = dataset.catalog[0].path
    with xr.open_dataset(path, engine="h5netcdf") as opened:
        if "x" in opened.coords:
            xv = opened["x"].values
            yv = opened["y"].values
        else:
            xv = opened["longitude"].values
            yv = opened["latitude"].values
    xs = slice(float(xv[0]), float(xv[-1]))
    ys = slice(float(yv[0]), float(yv[-1]))
    return xs, ys


def test_wind_interpolation_workflow():
    """Wind interpolation workflow using offline ``wind_3d_hourly_test`` fixtures (no CDS).

    Verifies:
    - Fixture dataset is registered and on disk
    - Model can be created and prepared
    - Capacity factor and wind-speed estimates run on the fixture extent
    """

    years = slice(2016, 2016)
    months = slice(1, 1)

    with Client(processes=True, threads_per_worker=1):
        ds_cls = load_dataset("wind_3d_hourly_test")
        ds = ds_cls(years=years, months=months)
        assert ds.downloaded, "Fixture NetCDF should be present"

        xs, ys = _fixture_xy_slices(ds)

        model = WindInterpolationModel(ds)
        assert model is not None
        model.prepare(force=True)

        turbine_name = "Enercon_E126_7500kW"

        cf_global = model.estimate(turbine=turbine_name)
        assert cf_global is not None
        assert isinstance(cf_global, (xr.DataArray, xr.Dataset))

        cf_region = model.estimate(turbine=turbine_name, xs=xs, ys=ys)
        assert cf_region is not None
        assert isinstance(cf_region, (xr.DataArray, xr.Dataset))

        speed = model.estimate(height=100.0, xs=xs, ys=ys)
        assert speed is not None
        assert isinstance(speed, xr.DataArray)

        cf_computed = cf_region.compute()
        assert cf_computed is not None
        max_cf = cf_computed.max()
        assert max_cf is not None
