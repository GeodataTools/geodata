# Copyright 2025 Michael Davidson (UCSD), Xiqiang Liu (UCSD), Keyu Long (UCSD)

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

"""Offline ERA5 datasets backed by committed NetCDF files under ``tests/fixtures/``.

On construction, small template files are **copied** into
``DATASET_ROOT_PATH / era5 / <weather_config> / …`` so paths stay compatible with
model code that uses :meth:`~geodata.model.results.BaseModelResult.ref_path`.

Importing this module registers ``wind_3d_hourly_test`` and ``wind_solar_hourly_test``
in :data:`geodata.datasets.registry`.
"""

from __future__ import annotations

import logging
import shutil
from pathlib import Path

from geodata.config import DATASET_ROOT_PATH

from .._base import AtomicDataset
from .wind_3d.hourly import ERA5Wind3DHourlyDataset
from .wind_solar.hourly import ERA5WindSolarHourlyDataset

logger = logging.getLogger(__name__)

# Paths must match tests/fixtures/era5/<config>/...
_FIXTURE_YEAR = 2016
_FIXTURE_MONTH = 1
_FIXTURE_DAY = 1


def _resolve_fixture_root(config_dirname: str) -> Path:
    """Return ``tests/fixtures/era5/<config_dirname>`` by walking parents of this file.

    Works for editable installs where the repo contains ``tests/fixtures``. Wheel-only
    installs without that tree raise ``FileNotFoundError``.
    """
    here = Path(__file__).resolve()
    for root in [here.parent, *here.parents]:
        candidate = root / "tests" / "fixtures" / "era5" / config_dirname
        if candidate.is_dir():
            return candidate
    raise FileNotFoundError(
        f"Could not find tests/fixtures/era5/{config_dirname} starting from {here}. "
        "Offline fixture datasets need the repository tests/fixtures tree (e.g. editable install)."
    )


def _copy_fixture_into_storage(template_root: Path, storage_root: Path, relative: Path) -> None:
    src = template_root / relative
    if not src.is_file():
        raise FileNotFoundError(f"Expected fixture NetCDF at {src}")
    dest = storage_root / relative
    dest.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dest)


class ERA5Wind3DHourlyTestDataset(ERA5Wind3DHourlyDataset):
    """Same schema as :class:`ERA5Wind3DHourlyDataset`, but points at a single local file.

    ``years`` / ``months`` passed to :meth:`__init__` do not expand the catalog; the
    catalog is always the fixture for ``{_FIXTURE_YEAR}/{_FIXTURE_MONTH:02d}/{_FIXTURE_DAY:02d}.nc``.
    """

    weather_config = "wind_3d_hourly_test"

    def _extra_setup(self, **kwargs):
        template_root = _resolve_fixture_root("wind_3d_hourly_test")
        self.storage_root = DATASET_ROOT_PATH / self.module / self.weather_config
        rel = (
            Path(str(_FIXTURE_YEAR))
            / f"{_FIXTURE_MONTH:02d}"
            / f"{_FIXTURE_DAY:02d}.nc"
        )
        _copy_fixture_into_storage(template_root, self.storage_root, rel)

    @property
    def catalog(self) -> list[AtomicDataset]:
        return [AtomicDataset(self, _FIXTURE_YEAR, _FIXTURE_MONTH, _FIXTURE_DAY)]

    def get_monthly_catalog(self, year: int, month: int) -> list[AtomicDataset]:
        """Only the committed fixture day exists under ``ref_path``; do not list full month."""
        if not isinstance(year, int):
            raise ValueError("year must be an integer")
        if not isinstance(month, int):
            raise ValueError("month must be an integer")
        if not 1 <= month <= 12:
            raise ValueError("month must be between 1 and 12")
        if not self.years.start <= year <= self.years.stop:
            raise ValueError(
                f"year must be between {self.years.start} and {self.years.stop}"
            )
        if not self.months.start <= month <= self.months.stop:
            raise ValueError(
                f"month must be between {self.months.start} and {self.months.stop}"
            )
        if year == _FIXTURE_YEAR and month == _FIXTURE_MONTH:
            return [AtomicDataset(self, year, month, _FIXTURE_DAY)]
        return []

    def _download_file(self, file: AtomicDataset):
        raise RuntimeError(
            f"{self.weather_config} uses committed fixtures under tests/fixtures; download is disabled."
        )


class ERA5WindSolarHourlyTestDataset(ERA5WindSolarHourlyDataset):
    """Same schema as :class:`ERA5WindSolarHourlyDataset`, but points at one monthly fixture file."""

    weather_config = "wind_solar_hourly_test"

    def _extra_setup(self, **kwargs):
        template_root = _resolve_fixture_root("wind_solar_hourly_test")
        self.storage_root = DATASET_ROOT_PATH / self.module / self.weather_config
        rel = Path(str(_FIXTURE_YEAR)) / f"{_FIXTURE_MONTH:02d}.nc"
        _copy_fixture_into_storage(template_root, self.storage_root, rel)

    @property
    def catalog(self) -> list[AtomicDataset]:
        return [AtomicDataset(self, _FIXTURE_YEAR, _FIXTURE_MONTH)]

    def _download_file(self, file: AtomicDataset):
        raise RuntimeError(
            f"{self.weather_config} uses committed fixtures under tests/fixtures; download is disabled."
        )


__all__ = [
    "ERA5Wind3DHourlyTestDataset",
    "ERA5WindSolarHourlyTestDataset",
]
