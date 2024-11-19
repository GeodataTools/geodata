# Supported Input/Output Formats

## Inputs

**Geodata** is currently optimized to work with the following datasets:

### ERA5

* [ERA5 hourly data on single levels](https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels?tab=overview)
* ERA5 3D level data [ADD 3D layer title and link]

### MERRA2

* [MERRA2 hourly, single-level surface flux diagnostics](https://disc.gsfc.nasa.gov/datasets/M2T1NXFLX_5.12.4/summary)
* [MERRA2 monthly mean, single-level surface flux diagnostics](https://disc.gsfc.nasa.gov/datasets/M2TMNXFLX_5.12.4/summary)
* [MERRA2 daily mean, single-level diagnostics](https://disc.gsfc.nasa.gov/datasets/M2SDNXSLV_5.12.4/summary)
* [MERRA2 hourly, single-level diagnostics](https://disc.gsfc.nasa.gov/datasets/M2T1NXSLV_5.12.4/summary)
* [MERRA2 hourly, single-level radiation diagnostics](https://disc.gsfc.nasa.gov/datasets/M2T1NXRAD_5.12.4/summary)
* [MERRA2 monthly mean, single-level diagnostics](https://disc.gsfc.nasa.gov/datasets/M2TMNXSLV_5.12.4/summary)
* [MERRA2 monthly mean, single-level radiation diagnostics](https://disc.gsfc.nasa.gov/datasets/M2TMNXRAD_5.12.4/summary)
* [MERRA2 hourly, single level aerosol diagnostics](https://disc.gsfc.nasa.gov/datasets/M2T1NXAER_5.12.4/summary)
* Functionality and documentation for additional datasets is planned.

Geodata can work with the following GIS files for masking:

* Tag Image File Format ([`.tif`](https://en.wikipedia.org/wiki/GeoTIFF))
* The shapefile format (`.shp`; geospatial vector data format for geographic information system software)

## Outputs

The following outputs are currently supported for climate data:

**Wind**

* Wind generation time-series ([MERRA2](../datasets/merra2/merra2_outputs.md#wind-generation-time-series), [ERA5](../datasets/era5/era5_outputs.md#wind-generation-time-series))
* Wind speed time-series ([MERRA2](../datasets/merra2/merra2_outputs.md#wind-speed-time-series), [ERA5](../datasets/era5/era5_outputs.md#wind-speed-time-series))
* Wind power density time-series ([MERRA2 only](../datasets/merra2/merra2_outputs.md#wind-power-density-time-series))


**Solar**

* Solar photovoltaic generation time-series ([MERRA2](../datasets/merra2/merra2_outputs.md#pv-generation-time-series), [ERA5](../datasets/era5/era5_outputs.md#solar-photovoltaic-generation-time-series))


**Temperature**

* Celsius Temperature (MERRA2 only) [ADD ERA5 temperature and links to both]


**Aerosols**

* PM2.5 time series ([MERRA2 only](../datasets/merra2/merra2_outputs.md#pm25-time-series))





