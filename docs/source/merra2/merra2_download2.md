
blah blah blah some test

Back to the [Table of Contents](https://github.com/east-winds/geodata/blob/master/doc/general/tableofcontents.md).

# Downloading MERRA2 Data and creating MERRA2 Cutouts

A short guide to use **geodata** to download MERRA2 data, and create cutouts - subsets of data based on specific time and geographic ranges - for [MERRA2 hourly, single-level surface flux diagnostics](https://disc.gsfc.nasa.gov/datasets/M2T1NXFLX_5.12.4/summary).

**Geodata** is currently optimized to work with the following MERRA2 datasets:
* [MERRA2 hourly, single-level surface flux diagnostics](https://disc.gsfc.nasa.gov/datasets/M2T1NXFLX_5.12.4/summary)
* [MERRA2 monthly mean, single-level surface flux diagnostics](https://disc.gsfc.nasa.gov/datasets/M2TMNXFLX_5.12.4/summary)
* [MERRA2 daily mean, single-level diagnostics](https://disc.gsfc.nasa.gov/datasets/M2SDNXSLV_5.12.4/summary)
* [MERRA2 hourly, single-level diagnostics](https://disc.gsfc.nasa.gov/datasets/M2T1NXSLV_5.12.4/summary)
* [MERRA2 hourly, single-level radiation diagnostics](https://disc.gsfc.nasa.gov/datasets/M2T1NXRAD_5.12.4/summary)
* [MERRA2 monthly mean, single-level diagnostics](https://disc.gsfc.nasa.gov/datasets/M2TMNXSLV_5.12.4/summary)
* [MERRA2 monthly mean, single-level radiation diagnostics](https://disc.gsfc.nasa.gov/datasets/M2TMNXRAD_5.12.4/summary)
* [MERRA2 hourly, single level aerosol diagnostics](https://disc.gsfc.nasa.gov/datasets/M2T1NXAER_5.12.4/summary)

## Download using the geodata package itself

Download methods for MERRA2 hourly data are built into the **geodata** package.  A basic example is as follows:

```
import logging
logging.basicConfig(level=logging.INFO)
import geodata

DS = geodata.Dataset(module="merra2",
					weather_data_config="surface_flux_hourly",
					years=slice(2012, 2012),
					months=slice(2,2))

if DS.prepared == False:
	DS.get_data()

DS.trim_variables(downloadedfiles = True)
```

Let's breakdown the code above:

```
import logging
logging.basicConfig(level=logging.INFO)
import geodata
```
Importing the logging package allows **geodata** to generate console input for download status, errors, etc.
Importing geodata is necessary to run **geodata**.

```
DS = geodata.Dataset(module="merra2",
					weather_data_config="surface_flux_monthly",
					years=slice(2012, 2012),
					months=slice(1,3))
```

`geodata.Dataset()` creates a dataset object via which you can download data files.  To create objects for MERRA2 hourly, single-level surface flux diagnostics, specify `module="merra2"`.  You must also have configured `merra2_dir` in `config.py` to point to a directory on your local machine  (to set up `config.py`, [see here](https://github.com/east-winds/geodata/blob/master/doc/general/packagesetup.md)).

The `years` and `months` parameters allow you to specify start years/months and end years/months for the data download.  The above example would download data for January to March 2012.  Ranges based on more granular time periods (such as day or hour) are not currently supported, but may be available in a future release.


Running the above code does not actually download the data yet.  Instead, it checks whether the indicated files for download are present in the local directory specified in `config.py`:

```
>> INFO:geodata.dataset:Directory /Users/johndoe/desktop/geodata/data/merra2 found, checking for completeness.
>> INFO:geodata.dataset:Directory complete.
```

and returns a `dataset` object indicating whether the data is "prepared."

```
<Dataset merra2 years=2012-2012 months=2-2 datadir=/Users/williamhonaker/desktop/davidson/data_for_geodata/merra2 Prepared>
```

A "prepared" dataset indicates that the directories for storing the data - which take the form `/merra2/{years}` (or `/merra2/{years}/{months}` for hourly data) for every unique time period in the data - have been created and are populated with downloaded data.  

```
if DS.prepared == False:
	DS.get_data()
```
If the dataset is "Unprepared", this conditional statement will download the actual netcdf files from GES DISC.

Files for the MERRA2 hourly, single-level surface flux diagnostics will be downloaded at the daily level, meaning that if you download data fo January 2011, 31 daily files will download to the directory `/merra2/2011/01`.  If you download data for January and February 2011, 31 daily files will download to the directory `/merra2/2011/01`, and 28 daily files will download to the directory `/merra2/2011/02`.

Each daily file is about 400mb in size, and contains all 46 variables detailed under the "Subset/Get Data" tab here: [MERRA2 hourly, single-level surface flux diagnostics](https://disc.gsfc.nasa.gov/datasets/M2T1NXFLX_5.12.4/summary).

Files for the MERRA2 monthly, single-level surface flux diagnostics will be downloaded at the monthly level, meaning that if you download data fo January 2011 to June 2011, 6 monthly files will download to the directory `/merra2/2011`.  These files contain the same variables as the hourly dataset described above.



```
DS.trim_variables(downloadedfiles = True)
```
To save hard disk space, **geodata** allows you to trim the downloaded datasets to just the variables needed for the package's supported outputs. Running the above code reduces file size from roughly 400mb per daily file to roughly 115mb per daily file by subsetting to only variables needed for generating wind and solar outputs.  

By default, `trim_variables()` subsets to both wind and solar data.  To further subset to just a single group of variables, you can pass `wind=False` or `solar=False` to exclude those groups of variables from the trimming process.

* Currently, only subsetting for MERRA2 wind data (`surface_flux_hourly`, `surface_flux_monthly`) is available.  Subsetting for solar data is planned for the near future.

For MERRA2 data, running `trim_variables()` as specified iterates over each downloaded file and subsets the data to the following wind variables indicated in the `surface_flux` configuration:

* ustar: surface_velocity_scale (m s-1)
* z0m: surface_roughness (m)
* disph: zero_plane_displacement_height (m)
* rhoa: air_density_at_surface (kg m-3)
* ulml: surface_eastward_wind (m s-1)
* vlml: surface_northward_wind (m s-1)
* tstar: surface_temperature_scale (K)
* hlml: surface_layer_height (m)
* tlml: surface_air_temperature (K)
* pblh: planetary_boundary_layer_height (m)
* hflux: sensible_heat_flux_from_turbulence (W m-2)
* eflux: total_latent_energy_flux (W m-2)


## Preparing the cutout

A cutout is the basis for any data or analysis output by the **geodata** package.  Cutouts are stored in the directory `cutout_dir` configured in `config.py` (to set up `config.py`, [see here](https://github.com/east-winds/geodata/blob/master/doc/general/packagesetup.md)).

```
cutout = geodata.Cutout(name="merra2-europe-sub24-2011-01",
                       module="merra2",
                       weather_data_config="surface_flux_hourly",
                       xs=slice(30, 41.56244222),
                       ys=slice(33.56459975, 35),
                       years=slice(2011, 2011),
                       months=slice(1,1))
```

To prepare a cutout, the following must be specified for `geodata.Cutout()`:

* The cutout name
* The source dataset
* Time range
* Geographic range as represented by `xy` coordinates.

The example in the code block above uses MERRA2 data, as specified by the `module` parameter.

```
module="merra2"
```

`xs` and `ys` in combination with the `slice()` function allow us to specify a geographic range based on longitude and latitude.  The above example subsets a portion of Europe.

`years` and `months` are used to subset the time range.  For both functions, the first value represents the start point, and the second value represents the end point.  The above example creates a cutout for January 2011.



## Creating a Cutout

To create a cutout, run:
```
cutout.prepare();
```
The **geodata** package will create a folder in the cutout directory (`cutout_dir`) you specified in `config.py` with the name specified in `geodata.Cutout()` (in the above example, `merra2-europe-sub24-2011-01`).  The folder, depending on the date range, will then contain one or more netcdf files containing data from the original MERRA2 files falling within the temporal and geographical ranges indicated when the cutout was created.

Before running `cutout.prepare()`, **geodata** will check for the presence of the cutout and abort if the cutout already exists.  If you want to force the regeneration of a cutout, run the command with the parameter `overwrite=True`.


## Cutout Metadata

You can query various metadata associated with a cutout.

```
cutout
```
returns the name, geographic and time range, and the preparation status of the cutout (i.e., whether cutout.prepare() has been run, creating the .nc files making up the cutout data).

```
<Cutout merra2-europe-sub24-2012-02 x=30.00-41.25 y=34.00-35.00 time=2012/2-2012/2 prepared>
```

`cutout.name` returns just the name:

```
'merra2-europe-sub24-2012-02'
```

`cutout.coords` returns coordinates:
```
Coordinates:
  * y           (y) float64 34.0 34.5 35.0
  * x           (x) float64 30.0 30.62 31.25 31.88 ... 39.38 40.0 40.62 41.25
  * time        (time) datetime64[ns] 2012-02-01T00:30:00 ... 2012-02-29T23:30:00
    lon         (x) float64 30.0 30.62 31.25 31.88 ... 39.38 40.0 40.62 41.25
    lat         (y) float64 34.0 34.5 35.0
  * year-month  (year-month) MultiIndex
  - year        (year-month) int64 2012
  - month       (year-month) int64 2
```

`cutout.meta` returns all associated metadata:

```
<xarray.Dataset>
Dimensions:     (time: 696, x: 19, y: 3, year-month: 1)
Coordinates:
  * y           (y) float64 34.0 34.5 35.0
  * x           (x) float64 30.0 30.62 31.25 31.88 ... 39.38 40.0 40.62 41.25
  * time        (time) datetime64[ns] 2012-02-01T00:30:00 ... 2012-02-29T23:30:00
    lon         (x) float64 30.0 30.62 31.25 31.88 ... 39.38 40.0 40.62 41.25
    lat         (y) float64 34.0 34.5 35.0
  * year-month  (year-month) MultiIndex
  - year        (year-month) int64 2012
  - month       (year-month) int64 2
Data variables:
    *empty*
Attributes:
    History:                           Original file generated: Wed Jul 23 15...
    Comment:                           GMAO filename: d5124_m2_jan10.tavg1_2d...
    Filename:                          MERRA2_400.tavg1_2d_flx_Nx.20120216.nc4
    Conventions:                       CF-1
    Institution:                       NASA Global Modeling and Assimilation ...
    References:                        http://gmao.gsfc.nasa.gov
    Format:                            NetCDF-4/HDF-5
    SpatialCoverage:                   global
    VersionID:                         5.12.4
    TemporalRange:                     1980-01-01 -> 2016-12-31
    identifier_product_doi_authority:  http://dx.doi.org/
    ShortName:                         M2T1NXFLX
    GranuleID:                         MERRA2_400.tavg1_2d_flx_Nx.20120216.nc4
    ProductionDateTime:                Original file generated: Wed Jul 23 15...
    LongName:                          MERRA2 tavg1_2d_flx_Nx: 2d,1-Hourly,Ti...
    Title:                             MERRA2 tavg1_2d_flx_Nx: 2d,1-Hourly,Ti...
    SouthernmostLatitude:              -90.0
    NorthernmostLatitude:              90.0
    WesternmostLongitude:              -180.0
    EasternmostLongitude:              179.375
    LatitudeResolution:                0.5
    LongitudeResolution:               0.625
    DataResolution:                    0.5 x 0.625
    Source:                            CVS tag: GEOSadas-5_12_4
    Contact:                           http://gmao.gsfc.nasa.gov
    identifier_product_doi:            10.5067/7MCPBJ41Y0K6
    RangeBeginningDate:                2012-02-16
    RangeBeginningTime:                00:00:00.000000
    RangeEndingDate:                   2012-02-16
    RangeEndingTime:                   23:59:59.000000
    module:                            merra2
```


To understand the variables that were downloaded, you can run:
```
cutout.dataset_module.weather_data_config
```

Which will return the selected download settings for the MERRA2 data.

```
{'surface_flux_hourly': {'tasks_func': <function geodata.datasets.merra2.tasks_monthly_merra2(xs, ys, yearmonths, prepare_func, **meta_attrs)>,
  'prepare_func': <function geodata.datasets.merra2.prepare_month_surface_flux(fn, year, month, xs, ys)>,
  'template': '/Users/johndoe/data_for_geodata/merra2/{year}/{month:0>2}/MERRA2_*.tavg1_2d_flx_Nx.*.nc4',
  'url': 'https://goldsmr4.gesdisc.eosdis.nasa.gov/data/MERRA2/M2T1NXFLX.5.12.4/{year}/{month:0>2}/MERRA2_{spinup}.tavg1_2d_flx_Nx.{year}{month:0>2}{day:0>2}.nc4',
  'fn': '/Users/johndoe/data_for_geodata/data_for_geodata/merra2/{year}/{month:0>2}/MERRA2_{spinup}.tavg1_2d_flx_Nx.{year}{month:0>2}{day:0>2}.nc4',
  'variables': ['ustar',
   'z0m',
   'disph',
   'rhoa',
   'ulml',
   'vlml',
   'tstar',
   'hlml',
   'tlml',
   'pblh',
   'hflux',
   'eflux']}}
```

To get just the variables that were downloaded for the `surface_flux_hourly` configuration, you can run:

```
cutout.dataset_module.weather_data_config['surface_flux_hourly']['variables']
```

Which will return a list of the variables included in the cutout when downloaded:

```
'ustar',
'z0m',
'disph',
'rhoa',
'ulml',
'vlml',
'tstar',
'hlml',
'tlml',
'pblh',
'hflux',
'eflux'
```
