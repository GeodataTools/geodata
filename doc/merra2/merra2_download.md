# Downloading MERRA2 Data for use with Geodata

A short guide to downloading data for MERRA2 data for use with **geodata**.

**Geodata** is currently optimized to work with the following MERRA2 datasets:
* [MERRA2 hourly, single-level surface flux diagnostics](https://disc.gsfc.nasa.gov/datasets/M2T1NXFLX_5.12.4/summary)
* [MERRA2 monthly mean, single-level surface flux diagnostics](https://disc.gsfc.nasa.gov/datasets/M2TMNXFLX_5.12.4/summary)

## Download using the geodata package itself

Download methods for MERRA2 hourly data are built into the **geodata** package.  A basic example is as follows:

```
import logging
logging.basicConfig(level=logging.INFO)
import atlite

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
import atlite
```
Importing the logging package allows **geodata** to generate console input for download status, errors, etc.
Importing atlite is necessary to run **geodata**.

```
DS = geodata.Dataset(module="merra2",
					weather_data_config="surface_flux_hourly",
					years=slice(2012, 2012),
					months=slice(2,2))
```

`geodata.Dataset()` creates a dataset object via which you can download data files.  To create objects for MERRA2 hourly, single-level surface flux diagnostics, specify `module="merra2"`.  You must also have configured `merra2_dir` in `config.py` to point to a directory on your local machine (see (insert link here) for details).

The `years` and `months` parameters allow you to specify start years/months and end years/months for the data download.  The above example would download data for February 2012.  Ranges based on more granular time periods (such as day or hour) are not currently supported, but may be available in a future release.


Running the above code does not actually download the data yet.  Instead, it checks whether the indicated files for download are present in the local directory specified in `config.py`:

```
>> INFO:geodata.dataset:Directory /Users/johndoe/desktop/geodata/data/merra2 found, checking for completeness.
>> INFO:geodata.dataset:Directory complete.
```

and returns a `dataset` object indicating whether the data is "prepared."

```
<Dataset merra2 years=2012-2012 months=2-2 datadir=/Users/williamhonaker/desktop/davidson/data_for_geodata/merra2 Prepared>
```

A "prepared" dataset indicates that the directories for storing the data - which take the form `/merra2/{years}/{months}` for every unique year-month combination in the data - have been created and are populated with downloaded data (Monthly-level data will download into just `/merra2/{years}`).  

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


To use downloaded data to create a cutout, see: [Creating Cutouts with MERRA2 Data (WIP)](https://github.com/east-winds/geodata/blob/master/doc/merra2_createcutout.md)
