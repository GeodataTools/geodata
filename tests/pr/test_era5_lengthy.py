# Copyright 2024 Xiqiang Liu

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

"""Tests in this file are lengthy due to the nature of the dataset being tested."""

import logging

import geodata

import xarray as xr

logging.basicConfig(level=logging.INFO)


# TODO: Test other functionalities with the 3D dataset
def get_data_configs() -> list[str]:
    return ["wind_3d_hourly"]


def get_bounds() -> list[list[int]]:
    return [[50, 0, 48, 3]]


def get_years() -> list[slice]:
    return [slice(2005, 2005)]


def get_months() -> list[slice]:
    return [slice(1, 2)]


def get_xs() -> list[slice]:
    return [slice(48.5, 49.5)]


def get_ys() -> list[slice]:
    return [slice(1, 2.5)]


def get_era5(data_config: str, bound: list[int], year: slice, month: slice):
    dataset = geodata.Dataset(
        module="era5",
        weather_data_config=data_config,
        years=year,
        months=month,
        bounds=bound,
    )
    if not dataset.prepared:
        dataset.get_data()
    return dataset


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

    for config, year, month, bound in zip(configs, years, months, bounds):
        dataset = get_era5(config, bound, year, month)
        dataset.trim_variables()
        for f in dataset.downloadedFiles:
            file_path = f[1]
            with xr.open_dataset(file_path) as ds:
                assert list(ds.data_vars) == dataset.weatherconfig["variables"]
