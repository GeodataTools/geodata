# Introduction to Geodata

The **geodata** package is a fork of [atlite](https://github.com/PyPSA/atlite), an xarray-based Python library for converting weather data (such as wind speeds, solar radiation, temperature and runoff) into power systems data (such as wind power, solar power, hydro power and heating demand time series). To learn more about **atlite** and its usage, see [the atlite readme](https://github.com/PyPSA/atlite/blob/master/README.rst) or its [Read the Docs page](https://atlite.readthedocs.io/en/latest/introduction.html).

# What is Geodata?

**Geodata** extracts global climate data from a given source and, for a user-specified time period and geographic range (subsets referred to as "cutouts"), converts it into a dataset of energy-related variables, such as wind generation or heat demand, for use in statistic analysis or other forms of academic research.

**Geodata** therefore provides an all-in-one Python interface for downloading, subsetting, and transforming climate data from large and cumbersome datasets into manageable data tables for use in regression, plotting, or other analyses.

# What inputs/outputs are currently possible with Geodata?

## Inputs
**Geodata** is currently optimized to work with the following datasets:
* [MERRA2 hourly, single-level surface flux diagnostics](https://disc.gsfc.nasa.gov/datasets/M2T1NXFLX_5.12.4/summary)
* [MERRA2 monthly mean, single-level surface flux diagnostics](https://disc.gsfc.nasa.gov/datasets/M2TMNXFLX_5.12.4/summary)
* [ERA5 hourly data on single levels](https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels?tab=overview)
* Functionality and documentation for additional datasets is planned.

## Outputs
The following outputs are currrently supported:
* Wind power generation


