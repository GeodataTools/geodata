# Copyright 2022 Xiqiang Liu

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
from itertools import product

import numpy as np
import xarray as xr

import geodata

logging.basicConfig(level=logging.INFO)


def get_data_configs() -> list[str]:
    return ["wind_solar_hourly"]


def get_bounds() -> list[list[int]]:
    return [[50, 0, 48, 3]]


def get_years() -> list[slice]:
    return [slice(2005, 2005)]


def get_months() -> list[slice]:
    return [slice(1, 1)]


def get_xs() -> list[slice]:
    return [slice(48.5, 49.5)]


def get_ys() -> list[slice]:
    return [slice(1, 2.5)]


def get_turbine() -> str:
    return "Suzlon_S82_1.5_MW"


def get_smooth() -> bool:
    return True


def get_panel() -> str:
    return "KANEKA"


def get_orientation() -> str:
    return "latitude_optimal"


def get_era5(data_config: str, bound: list[int], year: slice, month: slice):
    dataset = geodata.Dataset(
        module="era5",
        weather_data_config=data_config,
        years=year,
        months=month,
        bounds=bound,
    )
    if not dataset.prepared:
        dataset.get_data(testing=True)
    return dataset


def create_cutout(data_config: str, x: slice, y: slice, year: slice, month: slice):
    cutout = geodata.Cutout(
        name=f"{data_config}-europe-test-{year.start}-{month.start}",
        module="era5",
        weather_data_config=data_config,
        xs=x,
        ys=y,
        years=year,
        months=month,
    )
    cutout.prepare()
    return cutout


def create_wind_output(cutout: geodata.Cutout, turbine, smooth):
    ds_wind = geodata.convert.wind(cutout, turbine, smooth)
    df_wind = ds_wind.to_dataframe(name="wind")
    return df_wind


def create_windspd_output(cutout: geodata.Cutout, turbine):
    ds_windspd = geodata.convert.wind(cutout, turbine)
    df_windspd = ds_windspd.to_dataframe(name="windspd")
    return df_windspd


def create_pv_output(cutout: geodata.Cutout, panel, orientation):
    ds_pv = geodata.convert.pv(cutout, panel, orientation)
    df_pv = ds_pv.to_dataframe(name="pv")
    return df_pv


def test_download():
    configs = get_data_configs()
    years = get_years()
    months = get_months()
    bounds = get_bounds()

    for config, year, month, bound in zip(configs, years, months, bounds):
        dataset = get_era5(config, bound, year, month)
        assert dataset.prepared


def test_trim():
    configs = get_data_configs()
    years = get_years()
    months = get_months()
    bounds = get_bounds()

    for config, year, month, bound in product(configs, years, months, bounds):
        dataset = get_era5(config, bound, year, month)
        dataset.trim_variables()
        for f in dataset.downloadedFiles:
            file_path = f[1]
            with xr.open_dataset(file_path) as ds:
                assert list(ds.data_vars) == dataset.weatherconfig["variables"]


def test_cutout():
    configs = get_data_configs()
    years = get_years()
    months = get_months()
    xs = get_xs()
    ys = get_ys()

    for config, year, month, x, y in product(configs, years, months, xs, ys):
        cutout = create_cutout(config, x, y, year, month)
        assert cutout.prepared


def test_wind_output():
    configs = get_data_configs()
    years = get_years()
    months = get_months()
    xs = get_xs()
    ys = get_ys()

    supposed_dtypes = {
        "lat": np.floating,
        "lon": np.floating,
        "wind": np.floating,
    }

    for config, year, month, x, y in product(configs, years, months, xs, ys):
        cutout = create_cutout(config, x, y, year, month)
        turbine = get_turbine()
        smooth = get_smooth()
        df_wind = create_wind_output(cutout, turbine, smooth)

        for col in supposed_dtypes:
            assert np.issubdtype(df_wind[col].dtype, supposed_dtypes[col])


def test_windspd_output():
    configs = get_data_configs()
    years = get_years()
    months = get_months()
    xs = get_xs()
    ys = get_ys()

    supposed_dtypes = {
        "lat": np.floating,
        "lon": np.floating,
        "windspd": np.floating,
    }

    for config, year, month, x, y in product(configs, years, months, xs, ys):
        cutout = create_cutout(config, x, y, year, month)
        turbine = get_turbine()
        df_windspd = create_windspd_output(cutout, turbine)

        for col in supposed_dtypes:
            assert np.issubdtype(df_windspd[col].dtype, supposed_dtypes[col])


def test_pv_output():
    configs = [config for config in get_data_configs() if "hourly" in config]
    years = get_years()
    months = get_months()
    xs = get_xs()
    ys = get_ys()

    supposed_dtypes = {
        "lat": np.floating,
        "lon": np.floating,
        "pv": np.floating,
    }

    for config, year, month, x, y in zip(configs, years, months, xs, ys):
        cutout = create_cutout(config, x, y, year, month)
        panel = get_panel()
        orientation = get_orientation()
        df_pv = create_pv_output(cutout, panel, orientation)

        for col in supposed_dtypes:
            assert np.issubdtype(df_pv[col].dtype, supposed_dtypes[col])
