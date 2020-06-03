# Creating Cutouts with MERRA2 data

A short guide on how to use **geodata** to create cutouts - subsets of data based on specific time and geographic ranges - for [MERRA2 hourly, single-level surface flux diagnostics](https://disc.gsfc.nasa.gov/datasets/M2T1NXFLX_5.12.4/summary).

## Setup

To start, import the required dependencies:

```
import logging
logging.basicConfig(level=logging.INFO)
import atlite
```

`import atlite` is required to use **geodata**, while launching a logger allows for detailed debugging via the console.

## Preparing the cutout

A cutout is the basis for any data or analysis output by the **geodata** package.  Cutouts are stored in the directory `cutout_dir` configured in `config.py` (to set up `config.py`, see here LINK HERE)

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
{'surface_flux_hourly': {'tasks_func': <function atlite.datasets.merra2.tasks_monthly_merra2(xs, ys, yearmonths, prepare_func, **meta_attrs)>,
  'prepare_func': <function atlite.datasets.merra2.prepare_month_surface_flux(fn, year, month, xs, ys)>,
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