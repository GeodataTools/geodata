## Unreleased

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
- Added entry `surface_flux_dailymeans` for Daily Means Merra2 data (ADD LINK HERE).
- Added entries `slv_radiation_hourly`, `slv_radiation_monthly` for downloading and combining separate MERRA2 SLV and radition datasets into a single combined dataset for cutout create and solar output generation (ADD LINKS HERE).

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
- Added entry `surface_flux_monthly` for Monthly Merra2 data (ADD LINK HERE)
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

