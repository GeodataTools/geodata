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

import geodata

import xarray as xr
import numpy as np

logging.basicConfig(level=logging.INFO)


def get_data_configs() -> list[str]:
    return ["surface_flux_monthly"]


def get_bounds() -> list[list[int]]:
    return [[30, -10, 60, 10]]


def get_years() -> list[slice]:
    return [slice(2005, 2005)]


def get_months() -> list[slice]:
    return [slice(1, 3)]


def get_xs() -> list[slice]:
    return [slice(48.5, 49.5)]


def get_ys() -> list[slice]:
    return [slice(1, 2.5)]


def get_turbine() -> str:
    return "Suzlon_S82_1.5_MW"


def get_smooth() -> bool:
    return True


def get_var_height() -> str:
    return "lml"


def get_merra2(data_config: str, bound: list[int], year: slice, month: slice):
    dataset = geodata.Dataset(
        module="merra2",
        weather_data_config=data_config,
        years=year,
        months=month,
        bounds=bound,
    )
    if not dataset.prepared:
        dataset.get_data()
    return dataset


def create_cutout(data_config: str, x: slice, y: slice, year: slice, month: slice):
    cutout = geodata.Cutout(
        name="merra2-sample-01",
        module="merra2",
        weather_data_config=data_config,
        xs=x,
        ys=y,
        years=year,
        months=month,
    )
    cutout.prepare()
    return cutout


def create_wind_output(cutout: geodata.Cutout, turbine, smooth, var_height):
    ds_wind = geodata.convert.wind(cutout, turbine, smooth, var_height=var_height)
    df_wind = ds_wind.to_dataframe(name="wind")
    return df_wind


def create_windspd_output(cutout: geodata.Cutout, turbine, var_height):
    ds_windspd = geodata.convert.windspd(cutout, turbine=turbine, var_height=var_height)
    df_windspd = ds_windspd.to_dataframe(name="windspd")
    return df_windspd


def create_windwpd_output(cutout: geodata.Cutout, turbine, var_height):
    ds_windwpd = geodata.convert.windwpd(cutout, turbine=turbine, var_height=var_height)
    df_windwpd = ds_windwpd.to_dataframe(name="windwpd")
    return df_windwpd


def test_download():
    configs = get_data_configs()
    years = get_years()
    months = get_months()
    bounds = get_bounds()

    for config, year, month, bound in zip(configs, years, months, bounds):
        dataset = get_merra2(config, bound, year, month)
        assert dataset.prepared


def test_trim():
    configs = get_data_configs()
    years = get_years()
    months = get_months()
    bounds = get_bounds()

    for config, year, month, bound in zip(configs, years, months, bounds):
        dataset = get_merra2(config, bound, year, month)
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

    for config, year, month, x, y in zip(configs, years, months, xs, ys):
        cutout = create_cutout(config, x, y, year, month)
        assert cutout.prepared


def test_wind_output():
    configs = get_data_configs()
    years = get_years()
    months = get_months()
    xs = get_xs()
    ys = get_ys()

    for config, year, month, x, y in zip(configs, years, months, xs, ys):
        cutout = create_cutout(config, x, y, year, month)
        turbine = get_turbine()
        smooth = get_smooth()
        var_height = get_var_height()
        df_wind = create_wind_output(cutout, turbine, smooth, var_height)
        assert df_wind.dtypes.to_dict() == {
            "lon": np.dtype("float64"),
            "lat": np.dtype("float64"),
            "wind": np.dtype("float64"),
        }


def test_windspd_output():
    configs = get_data_configs()
    years = get_years()
    months = get_months()
    xs = get_xs()
    ys = get_ys()

    for config, year, month, x, y in zip(configs, years, months, xs, ys):
        cutout = create_cutout(config, x, y, year, month)
        turbine = get_turbine()
        var_height = get_var_height()
        df_windspd = create_windspd_output(cutout, turbine, var_height)
        assert df_windspd.dtypes.to_dict() == {
            "lon": np.dtype("float64"),
            "lat": np.dtype("float64"),
            "windspd": np.dtype("float32"),
        }


def test_windwpd_output():
    configs = get_data_configs()
    years = get_years()
    months = get_months()
    xs = get_xs()
    ys = get_ys()

    for config, year, month, x, y in zip(configs, years, months, xs, ys):
        cutout = create_cutout(config, x, y, year, month)
        turbine = get_turbine()
        var_height = get_var_height()
        df_windwpd = create_windwpd_output(cutout, turbine, var_height)
        assert df_windwpd.dtypes.to_dict() == {
            "lon": np.dtype("float64"),
            "lat": np.dtype("float64"),
            "windwpd": np.dtype("float32"),
        }
