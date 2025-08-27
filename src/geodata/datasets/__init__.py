# Copyright 2016-2017 Jonas Hoersch (FIAS), Tom Brown (FIAS), Markus Schlott (FIAS)
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
# along with this program. If not, see <http://www.gnu.org/licenses/>.


"""A interface to download and manage geospatial datasets.

Attributes:
    registry (dict): Registry of available datasets. You can retrieve relevant dataset classes using `registry.get(<weather_data_config>)`.
    DatasetType (type): A type to aid type hinting dataset classes.
"""

from . import era5, merra2
from ._base import DatasetType
from ._base import _registry as registry

__all__ = [
    "era5",
    "merra2",
    "registry",
    "register_hrrr",
    "DatasetType",
    "load_dataset",
    "list_datasets",
]


def load_dataset(weather_data_config: str) -> DatasetType:
    """Load a dataset based on the provided weather data configuration.
    This function retrieves the dataset class from the registry based on the provided configuration string.
    It is very similar to the `get` method of the registry (e.g. `registry.get(weather_data_config)`),
    but it is more explicit in its intent to load a dataset.

    Args:
        weather_data_config (str): The configuration string for the dataset.

    Returns:
        DatasetType: An instance of the dataset class corresponding to the configuration.

    Raises:
        ValueError: If the provided configuration is not registered in the registry.

    """
    if weather_data_config not in registry:
        raise ValueError(f"Dataset '{weather_data_config}' is not registered.")

    return registry[weather_data_config]


def list_datasets() -> list[str]:
    """List all registered datasets.

    Returns:
        list[str]: A list of dataset names that are registered in the registry.

    Raises:
        ValueError: If no datasets are registered in the registry.
    """

    if not registry:
        raise ValueError("No datasets are registered in the registry.")

    return list(registry.keys())


def register_hrrr():
    """Register the HRRR dataset with the registry. By default, the HRRR dataset is not registered to avoid import overhead with the HRRR Herbie package.
    This function can be called to register the HRRR dataset when needed."""
    from . import hrrr  # noqa: F401
