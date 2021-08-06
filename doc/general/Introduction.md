# Introduction to Geodata

**Geodata** is an `xarray`-enabled Python library of geospatial data collection and "pre-analysis" tools. Through the creation of shared scripts and documentation for analysis-ready physical variables, **geodata** streamlines the collection and use of geospatial datasets for natural science, engineering, and social science applications.

**Geodata** builds off the **[atlite](https://github.com/PyPSA/atlite)** library, which converts weather data (such as wind speeds, solar radiation, temperature and runoff) into power systems data (such as wind power, solar power, hydro power and heating demand time series). To learn more about atlite and its usage, see the atlite readme(https://github.com/PyPSA/atlite/blob/master/README.rst) or its documentation(https://atlite.readthedocs.io/en/latest/introduction.html).

# What is Geodata?

**Geodata** downloads global earth systems data from given sources and, for a user-specified time period and geographic range (subsets referred to as "cutouts"), converts it into a dataset of variables, such as wind generation or temperature.

**Geodata** therefore provides an all-in-one Python interface for downloading, subsetting, and transforming climate data from large and cumbersome datasets into manageable data tables for use in regression, plotting, or other analyses. Additionally, with a minimal amount of data consistency checks and metadata information, when one researcher goes through this exercise, everyone benefits.

# What inputs/outputs are currently possible with Geodata?

## Inputs

**Geodata** is currently optimized to work with the following datasets:

* [ERA5 hourly data on single levels](https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels?tab=overview)
* [MERRA2 hourly, single-level surface flux diagnostics](https://disc.gsfc.nasa.gov/datasets/M2T1NXFLX_5.12.4/summary)
* [MERRA2 monthly mean, single-level surface flux diagnostics](https://disc.gsfc.nasa.gov/datasets/M2TMNXFLX_5.12.4/summary)
* [MERRA2 daily mean, single-level diagnostics](https://disc.gsfc.nasa.gov/datasets/M2SDNXSLV_5.12.4/summary)
* [MERRA2 hourly, single-level diagnostics](https://disc.gsfc.nasa.gov/datasets/M2T1NXSLV_5.12.4/summary)
* [MERRA2 hourly, single-level radiation diagnostics](https://disc.gsfc.nasa.gov/datasets/M2T1NXRAD_5.12.4/summary)
* [MERRA2 monthly mean, single-level diagnostics](https://disc.gsfc.nasa.gov/datasets/M2TMNXSLV_5.12.4/summary)
* [MERRA2 monthly mean, single-level radiation diagnostics](https://disc.gsfc.nasa.gov/datasets/M2TMNXRAD_5.12.4/summary)
* [MERRA2 hourly, single level aerosol diagnostics](https://disc.gsfc.nasa.gov/datasets/M2T1NXAER_5.12.4/summary)
* Functionality and documentation for additional datasets is planned.

## Outputs

The following outputs are currrently supported:

### ERA5 Outputs:

**Wind**

* Wind generation time-series
* Wind speed time-series

**Solar**

* Solar photovoltaic generation time-series

### Merra2 Outputs:

**Wind**

* Wind generation time-series
* Wind speed time-series
* Wind power density time-series

**Solar**

* PV generation time-series

**Temperature**

* Celsius Temperature

**Aerosols**

* PM2.5 time series
