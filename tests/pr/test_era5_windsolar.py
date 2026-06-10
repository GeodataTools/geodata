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
from geodata.model.pvlib import Pvlib
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


def test_wind_solar_workflow():
    """Pvlib + wind interpolation using offline ``*_test`` fixtures (no CDS)."""

    years = slice(2016, 2016)
    months = slice(1, 1)

    with Client(processes=True, threads_per_worker=1):
        ds_cls = load_dataset("wind_solar_hourly_test")
        ds = ds_cls(years=years, months=months)
        assert ds.downloaded, "Fixture NetCDF should be present"

        xs, ys = _fixture_xy_slices(ds)

        model = Pvlib(ds)
        assert model is not None

        n_mods = 50
        n_strings = 1
        cec_modules = model.retrieve_sam("CECMod")
        module = cec_modules["Kaneka_U_SA105"]
        inv = model.retrieve_sam("CECInverter")["Fronius_USA__CL_33_3_Delta__208V_"]
        model.init_pv_system(
            arrays=None,
            surface_tilt=35,
            surface_azimuth=180,
            racking_model="open_rack",
            module_parameters=module,
            modules_per_string=n_mods,
            module_type="glass_polymer",
            module="Kaneka_U_SA105",
            strings_per_inverter=n_strings,
            inverter_parameters=inv,
        )
        assert model.pv_system is not None

        model.init_model_config(
            clearsky_model="haurwitz",
            transposition_model="perez",
            solar_position_method="nrel_numpy",
            airmass_model="kastenyoung1989",
            dc_model="cec",
            ac_model="sandia",
            aoi_model="physical",
            spectral_model="first_solar",
            dc_ohmic_model="no_loss",
        )
        assert model.config is not None

        ac_power_and_pv_capacity_global = model.estimate(
            years=years, months=months, xs=xs, ys=ys
        )
        assert ac_power_and_pv_capacity_global is not None
        assert isinstance(
            ac_power_and_pv_capacity_global, (xr.DataArray, xr.Dataset)
        )
        assert list(ac_power_and_pv_capacity_global.dims) == ["time", "x", "y"]

        wind_ds_cls = load_dataset("wind_3d_hourly_test")
        wind_ds = wind_ds_cls(years=years, months=months)
        assert wind_ds.downloaded, "Wind fixture NetCDF should be present"

        wxs, wys = _fixture_xy_slices(wind_ds)

        wind_model = WindInterpolationModel(wind_ds)
        wind_model.prepare()

        wind_speed = wind_model.estimate(
            years=years, months=months, xs=wxs, ys=wys, height=12
        )
        assert list(wind_speed.dims) == ["time", "x", "y"]
        assert (
            "valid_time" not in wind_speed.dims
            and "valid_time" not in wind_speed.coords
        )
