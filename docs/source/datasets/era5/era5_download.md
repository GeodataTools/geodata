# Download ERA5 Dataset and Create ERA5 Cutouts

A short guide to using **geodata** to downloading and creating cutouts from ERA5 data from the [Copernicus Data Store](https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels?tab=overview).

## Download Dataset

Download methods for ERA5 data are built into the **geodata** package.  A basic example is as follows:


To start, import the required dependencies:

```python
import geodata

dataset = geodata.Dataset(
    module="era5",
    weather_data_config="wind_solar_hourly",
    years=slice(2005, 2005),
    months=slice(1, 2),
    bounds=[50, -3, 45, 3],
)

if not dataset.prepared:
    dataset.get_data()

dataset.trim_variables(downloadedfiles=True)
```

Let's breakdown the code above:

```python
import geodata
```
Importing the logging package allows **geodata** to generate console input for download status, errors, etc.
Importing geodata is necessary to run **geodata**.

```python
dataset = geodata.Dataset(
    module="era5",
    weather_data_config="wind_solar_hourly",
    years=slice(2005, 2005),
    months=slice(1, 2),
    bounds=[50, -3, 45, 3],
)
```

`geodata.Dataset()` creates a dataset object via which you can download data files.  To create objects for ERA5 data, specify `module="era5"`.  You must also have configured `era5_dir` in `config.py` to point to a directory on your local machine  (to set up `config.py`, [see here](../../quick_start/packagesetup.md)).

The `years` and `months` parameters allow you to specify start years/months and end years/months for the data download.  The above example would download data for January to February 2005.  Ranges based on more granular time periods (such as day or hour) are not currently supported, but may be available in a future release.

Running the above code does not actually download the data yet.  Instead, it checks whether the indicated files for download are present in the local directory specified in `config.py`:

```
>> INFO:geodata.dataset:Directory /home/user/.local/geodata/era5 found, checking for completeness.
>> INFO:geodata.dataset:Directory complete.
```

and returns a `dataset` object indicating whether the data is "prepared."

```
<Dataset era5 years=2005-2005 months=1-2 datadir=/home/user/.local/geodata/era5 Prepared>
```

A "prepared" dataset indicates that the directories for storing the data - which take the form `/era5/{years}/{months}` for every unique time period in the data - have been created and are populated with downloaded data.  

```python
if not dataset.prepared:
    dataset.get_data()
```
If the dataset is "Unprepared", this conditional statement will download the actual netcdf files from the ERA5 CDS API.


```python
dataset.trim_variables(downloadedfiles = True)
```
To save hard disk space, **geodata** allows you to trim the downloaded datasets to just the variables needed for the package's supported outputs.

## Create Cutout

A cutout is the basis for any data or analysis output by the **geodata** package.  Cutouts are stored in the directory `cutout_dir` configured in `config.py` (to set up `config.py`, [see here](../../quick_start/packagesetup.md)).

```python
cutout = geodata.Cutout(
    name="era5-europe-test-2011-01",
    module="era5",
    weather_data_config="era5_monthly",
    xs=slice(30, 41.56244222),
    ys=slice(33.56459975, 35),
    years=slice(2011, 2011),
    months=slice(1, 1),
)
```

To prepare a cutout, the following must be specified for `geodata.Cutout()`:

* The cutout name
* The source dataset
* Time range 
* Geographic range as represented by `x` and `y` coordinates.

The example in the code block above uses ERA5 data, as specified by the `module=era5` parameter.

`xs` and `ys` in combination with the `slice()` function allow us to specify a geographic range based on longitude and latitude.  The above example subsets a portion of Europe.

`years` and `months` are used to subset the time range.  For both functions, the first value represents the start point, and the second value represents the end point.  The above example creates a cutout for January 2011.

## Prepare Cutout

To create a cutout, we can run `cutout.prepare()`. The **geodata** package will create a folder in the cutout directory you specified in `config.py` with the name specified in `geodata.Cutout()` (in the above example, `era5-europe-test-2011-01`).  The folder, depending on the date range, will then contain one or more monthly netcdf files containing ERA5 data corresponding to the temporal and geographical ranges indicated when the cutout was created.  Data files in the cutout folder will be at the monthly level - i.e., there will be one file for each month in the specified download time range.

To actually create the cutout, you must run `cutout.prepare()`.  Upon running `cutout.prepare()`, **geodata** will check for the presence of the cutout and abort if the cutout already exists.  If you want to force the regeneration of a cutout, run the command with the parameter `overwrite=True`.


## Cutout Metadata

You can query various metadata associated with a cutout. Querying returns the name, geographic and time range, and the preparation status of the cutout (i.e., whether cutout.prepare() has been run, creating the .nc files making up the cutout data).

```python
<Cutout era5-europe-test-2011-01 x=30.00-41.50 y=34.81-33.81 time=2011/1-2011/1 prepared>
```

- `cutout.name` returns just the name:

```
'era5-europe-test-2011-01'
```

- `cutout.coords` returns coordinates:
```python
Coordinates:
  * x           (x) float32 30.0 30.25 30.5 30.75 31.0 ... 40.75 41.0 41.25 41.5
  * y           (y) float32 34.815 34.565 34.315 34.065 33.815
  * time        (time) datetime64[ns] 2011-01-01 ... 2011-01-31T23:00:00
    lon         (x) float32 ...
    lat         (y) float32 ...
  * year-month  (year-month) MultiIndex
  - year        (year-month) int64 2011
  - month       (year-month) int64 1
```

- `cutout.meta` returns all associated metadata:
```python
<xarray.Dataset>
Dimensions:     (time: 744, x: 47, y: 5, year-month: 1)
Coordinates:
  * x           (x) float32 30.0 30.25 30.5 30.75 31.0 ... 40.75 41.0 41.25 41.5
  * y           (y) float32 34.815 34.565 34.315 34.065 33.815
  * time        (time) datetime64[ns] 2011-01-01 ... 2011-01-31T23:00:00
    lon         (x) float32 ...
    lat         (y) float32 ...
  * year-month  (year-month) MultiIndex
  - year        (year-month) int64 2011
  - month       (year-month) int64 1
Data variables:
    height      (y, x) float32 ...
Attributes:
    Conventions:  CF-1.6
    history:      2020-03-14 04:29:35 GMT by grib_to_netcdf-2.16.0: /opt/ecmw...
    module:       era5
    view:         {'x': slice(30, 41.56244222, None), 'y': slice(35, 33.56459...
```

To understand the variables that were downloaded, you can run:
```
cutout.dataset_module.weather_data_config
```

Which will return the selected download settings for the ERA5 data.
<!-- We need to update this after ERA5 keyword/variable renaming -->

```python
{'wind_solar_hourly': {'api_func': <function geodata.datasets.era5.api_hourly_era5(toDownload, bounds, download_vars, product, product_type)>,
  'file_granularity': 'monthly',
  'tasks_func': <function geodata.datasets.era5.tasks_monthly_era5(xs, ys, yearmonths, prepare_func, **meta_attrs)>,
  'meta_prepare_func': <function geodata.datasets.era5.prepare_meta_era5(xs, ys, year, month, template, module, **kwargs)>,
  'prepare_func': <function geodata.datasets.era5.prepare_month_era5(fn, year, month, xs, ys)>,
  'template': '/Users/johndoe/data_for_geodata/era5/{year}/{month:0>2}/wind_solar_hourly.nc',
  'fn': '/Users/johndoe/data_for_geodata/era5/{year}/{month:0>2}/wind_solar_hourly.nc',
  'product': 'reanalysis-era5-single-levels',
  'product_type': 'reanalysis',
  'variables': ['100m_u_component_of_wind',
   '100m_v_component_of_wind',
   '2m_temperature',
   'runoff',
   'soil_temperature_level_4',
   'surface_net_solar_radiation',
   'surface_pressure',
   'surface_solar_radiation_downwards',
   'toa_incident_solar_radiation',
   'total_sky_direct_solar_radiation_at_surface',
   'forecast_surface_roughness',
   'orography'],
  'meta_attrs': {'Conventions': 'CF-1.6',
   'history': '2022-01-11 00:43:10 GMT by grib_to_netcdf-2.23.0: /opt/ecmwf/mars-client/bin/grib_to_netcdf -S param -o /cache/data3/adaptor.mars.internal-1641861772.7584445-3425-17-8ff8bcdc-2c08-4aed-95db-e0337692a384.nc /cache/tmp/8ff8bcdc-2c08-4aed-95db-e0337692a384-adaptor.mars.internal-1641861207.3577378-3425-30-tmp.grib',
   'module': 'era5'}},
 'wind_solar_monthly': {'api_func': <function geodata.datasets.era5.api_monthly_era5(toDownload, bounds, download_vars, product, product_type)>,
  'file_granularity': 'monthly',
  'tasks_func': <function geodata.datasets.era5.tasks_monthly_era5(xs, ys, yearmonths, prepare_func, **meta_attrs)>,
  'meta_prepare_func': <function geodata.datasets.era5.prepare_meta_era5(xs, ys, year, month, template, module, **kwargs)>,
  'prepare_func': <function geodata.datasets.era5.prepare_month_era5(fn, year, month, xs, ys)>,
  'template': '/Users/johndoe/data_for_geodata/era5/{year}/{month:0>2}/wind_solar_monthly.nc',
  'fn': '/Users/johndoe/data_for_geodata/era5/{year}/{month:0>2}/wind_solar_monthly.nc',
  'product': 'reanalysis-era5-single-levels-monthly-means',
  'product_type': 'monthly_averaged_reanalysis',
  'variables': ['100m_u_component_of_wind',
   '100m_v_component_of_wind',
   '2m_temperature',
   'runoff',
   'soil_temperature_level_4',
   'surface_net_solar_radiation',
   'surface_pressure',
   'surface_solar_radiation_downwards',
   'toa_incident_solar_radiation',
   'total_sky_direct_solar_radiation_at_surface',
   'forecast_surface_roughness',
   'orography']}}
```

