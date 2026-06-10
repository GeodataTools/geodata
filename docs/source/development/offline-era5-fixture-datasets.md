# Offline ERA5 fixture datasets (`*_test` weather configs)

This document records the **design and implementation plan** for small, committed NetCDF fixtures used in automated tests—without calling the CDS API or relying on `DATASET_ROOT_PATH` downloads.

## Goals

- Ship **minimal** ERA5-shaped files in the repository for CI and local testing.
- Expose them via **`load_dataset("…_test")`** so code paths mirror production (`wind_3d_hourly`, `wind_solar_hourly`) while staying **offline**.
- Avoid coupling tests to arbitrary year/month ranges: fixture datasets should use a **fixed catalog** (typically a single file) even if `BaseDataset.__init__` still requires `years` / `months` arguments (those values can be **ignored** for catalog construction in test configs).

## Non-goals

- The legacy **`geodata.dataset.Dataset`** (`module=` + `weather_data_config=` dict) is **not** in scope; the plan targets **`load_dataset` + `BaseDataset` subclasses** used by models and current tests.

## Current fixture layout (repository)

Fixtures live under **`tests/fixtures/`** so they stay close to pytest and do not inflate the installable package unless explicitly packaged later.

| Test weather config (planned) | Mirrors production config | On-disk layout under `tests/fixtures/era5/` |
|-------------------------------|---------------------------|---------------------------------------------|
| `wind_3d_hourly_test` | `wind_3d_hourly` (`frequency="daily"`) | `wind_3d_hourly_test/2016/01/01.nc` |
| `wind_solar_hourly_test` | `wind_solar_hourly` (default `frequency="monthly"`) | `wind_solar_hourly_test/2016/01.nc` |

Production datasets store files under:

`DATASET_ROOT_PATH / <module> / <weather_config> / …`

with:

- **Daily** (3D wind): `…/<year>/<month>/<day>.nc`
- **Monthly** (wind/solar hourly): `…/<year>/<month>.nc`

The fixture tree **matches those relative paths** so `AtomicDataset.path` resolution stays aligned with the real datasets.

## Registry and naming

- Each test variant is a **`BaseDataset` subclass** with `weather_config = "wind_3d_hourly_test"` or `"wind_solar_hourly_test"`.
- Subclasses are registered automatically via `BaseDataset.__init_subclass__` into `geodata.datasets.registry`.
- Callers use **`load_dataset("wind_3d_hourly_test")`** (same pattern as production).

## Behavioral contract

### Storage root

Fixture classes should set **`storage_root`** to the directory that contains the fixture tree for that config—for example, the absolute path to `tests/fixtures/era5/wind_3d_hourly_test` resolved at runtime (repo-relative or via `importlib.resources` if fixtures are ever packaged).

### Catalog

Override **`catalog`** so it returns **only** the `AtomicDataset` entries that refer to committed files (commonly **one** file):

- 3D wind: one daily file, e.g. `(year=2016, month=1, day=1)` → `…/2016/01/01.nc`
- Wind/solar: one monthly file, e.g. `(year=2016, month=1)` → `…/2016/01.nc`

Constructor arguments **`years` / `months`** may remain required by `BaseDataset.__init__` but **need not drive** the fixture catalog.

### Download

- **`download()`** must **not** call CDS: implement as a no-op or raise a clear error if invoked.
- **`_download_file`** should not perform network I/O.

### Prepared state

`downloaded` should become **`True`** when fixture files exist (the default `_check_downloaded()` loop over `catalog` is sufficient if paths resolve correctly).

## Models and `SUPPORTED_WEATHER_DATA_CONFIGS`

`BaseModel` validates both **`weather_config`** and **`source.downloaded`**. Any model that should run on fixtures must **allow** the `*_test` config names—e.g. extend `SUPPORTED_WEATHER_DATA_CONFIGS` on `WindInterpolationModel`, pvlib-related models, and any other entry points used in tests—to include `wind_3d_hourly_test` / `wind_solar_hourly_test` (or document a single shared alias strategy).

## Implementation checklist

1. Add **`ERA5Wind3DHourlyTestDataset`** / **`ERA5WindSolarHourlyTestDataset`** (names may vary) beside the existing ERA5 hourly classes, or in a small `fixture.py` module imported from `era5` packages so subclasses register on import.
2. Wire **`storage_root`** to `tests/fixtures/era5/<config>/` (resolve path robustly from the repo root or test layout).
3. Override **`catalog`** to the fixed fixture file(s); ignore user `years`/`months` for catalog purposes (documented).
4. Override **`download`** / **`_download_file`** to prevent CDS usage.
5. Update **`SUPPORTED_WEATHER_DATA_CONFIGS`** on affected models.
6. Add or adjust tests: `load_dataset("…_test")`, assert `downloaded`, **no** `download()`, then run the intended model or pipeline assertion.

## References (code)

- Registry: `geodata.datasets._base.BaseDataset.__init_subclass__`
- Paths: `AtomicDataset.path` in `geodata.datasets._base`
- Legacy downloader: `geodata.dataset.Dataset` (separate from this plan)
