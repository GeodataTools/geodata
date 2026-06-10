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

# import logging

# from geodata.datasets import DatasetType, load_dataset

# logging.basicConfig(level=logging.INFO)


# def get_data_configs() -> list[str]:
#     return [
#         "surface_flux_monthly",
#         "slv_radiation_monthly",
#         "surface_flux_hourly",
#         "slv_radiation_hourly",
#     ]


# def get_bounds() -> list[list[int]]:
#     return [[30, -10, 60, 10]]


# def get_years() -> list[slice]:
#     return [slice(2005, 2005)]


# def get_months() -> list[slice]:
#     return [slice(1, 1)]


# def get_merra2(data_config: str, bound: list[int], year: slice, month: slice):
#     dataset_cls = load_dataset(data_config)
#     dataset: DatasetType = dataset_cls(
#         years=year, months=month, bounds=bound, testing=True
#     )
#     if not dataset.downloaded:
#         dataset.download()
#     return dataset


# def test_download():
#     configs = get_data_configs()
#     years = get_years()
#     months = get_months()
#     bounds = get_bounds()

#     for config, year, month, bound in zip(configs, years, months, bounds):
#         dataset = get_merra2(config, bound, year, month)
#         assert dataset.downloaded
