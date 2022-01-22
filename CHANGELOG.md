## Unreleased

## datasets/
* Removal of `sarah.py`

## geodata/config-default.py
* Removal of `sarah_dir`

## era5.py
* Remove following unused functions:
  - get_data
  - _add_height
  - _add_area
  - prepare_for_sarah
* Remove various unused import statements

## dataset.py
### get_data()
* Add `testing` flag to enable download of only first entry in dataset download list for debugging purposes.

## preparation.py
* Added lines to close tempfiles.  This was causing Win32 permissions issues for windows machines.

## doc/
* Various updates to setup documentation regarding **rasterio** installation, Merra2 API access, and Windows directory formats.

## merra2.py
* Added lines to close tempfiles used in `api_merra2` function.  This was causing Win32 permissions issues for windows machines.
* Move function lines following API request inside `try` block so that "Success!" info message does not fire in case of HTTP error.

## geodata.egg-info/
* Remove `geodata.egg-info/` from repo.  This should prevent the folder from being tracked.  Otherwise the presence of the folder in the master report trumps listing the file in `.gitignore`.

## preparation.py
* Fixed issue where `meta.attrs.setdefault()` could attempt to assign value to subset of string in case of meta file being read in from `xarray.opendataset()`.  Safer to redeclare `view` key as property and then use `meta.attrs.setdefault()`.
https://github.com/east-winds/geodata/pull/18

## Merra2 PM25 Update
PR: https://github.com/east-winds/geodata/pull/14

### tests
* Added `merra_aerosol_hourly_test.ipynb` to test aerosol hourly config for MERRA2.
* Removed locally-generated output text from jupyter notebook files in `tests` folder.

### merra2.py
* [`weather_data_config`]
- Added entry `surface_aerosol_hourly` for [MERRA2 hourly, single level aerosol diagnostics](https://disc.gsfc.nasa.gov/datasets/M2T1NXAER_5.12.4/summary)

### dataset.py
* Removed extra text 'pi' from "files not found" message (lines 117, 137)

### convert.py
* Added `pm25` and `convert_pm25` for generating pm25 time series output.

### geodata.egg-info/PKG-INFO
* Updated package info.

## Merra2 Daily Update

### tests
* Renamed existing MERRA2 tests to `merra2_surface_flux_hourly_test.ipynb` and `merra2_surface_flux_monthly_test.ipynb`.
* Added `merra2_surface_flux_dailymeans_test.ipynb` to test daily means config for MERRA2.
* Added `merra2_slv_radiation_hourly_test.ipynb` and `merra2_slv_radiation_monthly_test.ipynb` to test new slv and radiatian MERRA2 config.
* Added `era5_test.ipynb`.

### dataset.py
* [`Dataset`]
- Added case for `file_granularity=='dailymeans'` to initialization.
- Added cases `file_granularity=='daily_multiple'` and `file_granularity=='monthly_multiple'` to allow download from both MERRA2 SLV and radiation datasets and combination into a single combined dataset.

* [`get_data`]
- Added case for `file_granularity=='daily_multiple'` and `file_granularity=='monthly_multiple'`

### preparation.py
* [`cutout_get_meta`] Added in timestamp creation case for `daily means` MERRA2 data.

### merra2.py
* [`weather_data_config`]
- Added entry `surface_flux_dailymeans` for Daily Means Merra2 data - [MERRA2 daily mean, single-level diagnostics](https://disc.gsfc.nasa.gov/datasets/M2SDNXSLV_5.12.4/summary).
- Added entries `slv_radiation_hourly`, `slv_radiation_monthly` for downloading and combining separate MERRA2 SLV and radition datasets into a single combined dataset for cutout create and solar output generation:
    * `slv_radiation_hourly`:
        - [MERRA2 hourly, single-level diagnostics](https://disc.gsfc.nasa.gov/datasets/M2SDNXSLV_5.12.4/summary)
        - [MERRA2 hourly, single-level radiation diagnostics](https://disc.gsfc.nasa.gov/datasets/M2T1NXRAD_5.12.4/summary)
    * `slv_radiation_monthly`:
        - [MERRA2 monthly mean, single-level diagnostics](https://disc.gsfc.nasa.gov/datasets/M2TMNXSLV_5.12.4/summary)
        - [MERRA2 monthly mean, single-level radiation diagnostics](https://disc.gsfc.nasa.gov/datasets/M2TMNXRAD_5.12.4/summary)

* Other functions
- Added prepare function `prepare_dailymeans_surface_flux`.
- Added prepare function `prepare_slv_radiation`.


## Merra2 Monthly Update

### tests
* Added `merra2_hourly_test.ipynb` and `merra2_monthly_test.ipynb` for simple functionality testing.

### cutout.py
* [`cutoutparams`]
- Added `weather_data_config` to `cutoutparams` in order to allow per-config specification of cutout file granularity and prepare function. `weather_data_config` for a given cutout should be the same as the dataset used to download the source data files.  A future improvement could make this config parameter auto-populate from the source datafiles.
- Updated `meta_data_config` property so that it is built directly from the module's `weather_data_config`.
- Added `file_granularity` to `meta_data_config` property so that proper file granularity can be accounted for when defining tasks for cutout preparation.

### dataset.py
* [`Dataset`]
- Added param `weather_data_config` to object initialization.  Dataset objects now to be associated with a given `config` as defined in the `weather_data_confg` of the relevant data source module. Config name accessed at `self.config`, config details accessed at `self.weatherconfig`.
- Added param `totalFiles` to keep list of all data source files associated with Dataset object, regardless of `self.prepared` value.
- When intializing list of data source filenames, Dataset now looks for `file_granularity` property in specified `weather_data_config`, rather than static value `dataset_module.daily_files`.
- Cleaned up for loops for defining filenames to downlad (lines 85-118).  Both loops now work on a per-`config` basis rather than trying to account for multiple datasource in a single object.  Both loops now also only iterate through year-month tuples (monthly config) or year-month-day tuples (for daily) instead of nested loops through extra `config` setings.
- Cleaned up log output in loops (lines 96, 112) to print only name of each missing file (previously was printing file name once per number of properties in `weather_data_confg`).
- Removed references to `savedFiles` in favor of `downloadedFiles`.

* [`get_data()`]
- Testing parameter removed in favor of future functionality around daily downloads (ie, no need for an explicit option to download only one file).
- Download for loop now only uses filename list (`self.toDownload`) as defined on initialization of Dataset object - ie now only downloads files from a single data source as defined by the config defined in `self.weatherconfig`. Method previously tried to account for files across multiple sources (ie, a wind source and a solar source).
- `self.prepared` now auto-updates to `true` when `get_data()` successfully downloads all files.

* [`trim_variables()`]
- Now only refers to `downloadedFiles`.
- Now only refers to variables defined in Dataset's `weatherconfig`.

### merra2.py
* [`weather_data_config`]
- Added entry `surface_flux_monthly` for Monthly Merra2 data - [MERRA2 monthly mean, single-level surface flux diagnostics](https://disc.gsfc.nasa.gov/datasets/M2TMNXFLX_5.12.4/summary)
- To all entry added `file_granularity` property to indicated time granularity for file download.
* Preparation functions
- Separated cutout task preparation function into `tasks_monthly_merra2` and `tasks_daily_merra2` to account for differences in file structure between daily and hourly file granularities.

### preparation.py
* [123-131] - List of tasks for cutout preparation now accounts for differences between file granularites as specifiec in source `weather_data_confg`.
* [218-233] - Timestamp construction now depending on source config `file_granularity`.


## Documentation - 04/10/20

### ERA5
* [`era5_download.md`] Add guide for downloading and creating cutouts from ERA5 data.
* [`era5_example_output.md`] Add example guide for generating output for ERA5 data.
* [`era5_outputs.md`] Add guide for outputs for ERA5.
* [`era5.ipynb`] Add jupyter notebook to walk through entire ERA5 download, cutout, and output process.

### MERRA2
* [`merra2_outputs.md`] Roll up merra2 output guide files into a single file.


## Documentation - 04/04/20
* [`merra2_example_output.md`] Add example guide for generating output for MERRA2.
* [`merra2_configs_and_outputs.md`] Add guide for configs and outputs for MERRA2.
* [`merra2_output_wind.md`] Add guide for wind outputs in MERRA2.
* [`merra2_output_solar.md`] Add starter guide for solar outputs in MERRA2.
* [`merra2.ipynb`] Add jupyter notebook to walk through entire MERRA2 download, cutout, and output process.


## Documentation - 03/28/20
* [`Introduction.md`] Add introduction to the package.
* [`tableofcontents.md`] Add list of links to relevant documentation.
* [`packagesetup.md`] Add guide for installing package.
* [`merra2_setup.md`] Add guide for setting up access to merra2 hourly surface flux data.
* [`merra2_download.md`] Add guide for downloading merra2 data.
* [`merra2_createcutout.md`] Add guide for creating cutouts for merra2 data.
* [`era5_setup.md`] Add guide for setting up access to ERA5 data from CDS.

## Wind parameterizations - 11/29/19
* [`wind.py`] Add wind parameterizations (given by `extrap_fn`) for different stability correction functions
* [`doc/parameterizations/wind_parameterizations.tex.md`] Add citations and formulae for wind parameterizations (LaTeX)
* [`doc/parameterizations/wind_parameterizations.md`] Citations and formulae for wind parameterizations (created automatically from `wind_parameterizations.tex.md`)

