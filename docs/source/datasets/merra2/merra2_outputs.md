# MERRA2 Ouputs Documentation

This notebook includes these following sections:

- [MERRA2 Configs and Outputs](https://github.com/east-winds/geodata/blob/master/doc/merra2/merra2_outputs.md#merra2-configs-and-outputs)
- [Generating Wind Outputs with MERRA2](https://github.com/east-winds/geodata/blob/master/doc/merra2/merra2_outputs.md#generating-wind-outputs-with-merra2-data)
- [Generating Solar Outputs with MERRA2](https://github.com/east-winds/geodata/blob/master/doc/merra2/merra2_outputs.md#generating-solar-outputs-with-merra2-data)
- [Generating Aerosol Outputs with MERRA2](https://github.com/east-winds/geodata/blob/master/doc/merra2/merra2_outputs.md#generating-aerosol-outputs-with-merra2-data)
- [An Example Output](https://github.com/east-winds/geodata/blob/master/doc/merra2/merra2_outputs.md#example-output)


# MERRA2 Configs and Outputs

In this section, we provide a list of currently possible **geodata** configs and outputs for MERRA2 data from NASA's GES DISC.

## Dataset Configurations

**geodata** is currently optimized to work with the following MERRA2 datasets:

* `surface_flux_hourly`: [MERRA2 hourly, single-level surface flux diagnostics](https://disc.gsfc.nasa.gov/datasets/M2T1NXFLX_5.12.4/summary)
* `surface_flux_monthly`: [MERRA2 monthly mean, single-level surface flux diagnostics](https://disc.gsfc.nasa.gov/datasets/M2TMNXFLX_5.12.4/summary)
* `surface_flux_dailymeans`: [MERRA2 daily mean, single-level diagnostics](https://disc.gsfc.nasa.gov/datasets/M2SDNXSLV_5.12.4/summary)
* `slv_radiation_hourly`: 
    - [MERRA2 hourly, single-level diagnostics](https://disc.gsfc.nasa.gov/datasets/M2T1NXSLV_5.12.4/summary)
    - [MERRA2 hourly, single-level radiation diagnostics](https://disc.gsfc.nasa.gov/datasets/M2T1NXRAD_5.12.4/summary)
* `slv_radiation_monthly`:
    - [MERRA2 monthly mean, single-level diagnostics](https://disc.gsfc.nasa.gov/datasets/M2TMNXSLV_5.12.4/summary)
    - [MERRA2 monthly mean, single-level radiation diagnostics](https://disc.gsfc.nasa.gov/datasets/M2TMNXRAD_5.12.4/summary)
* `surface_aerosol_hourly`: [MERRA2 hourly, single level aerosol diagnostics](https://disc.gsfc.nasa.gov/datasets/M2T1NXAER_5.12.4/summary)


## Variable Configurations

MERRA2 surface flux data currently supports the following variable configs:

**Wind** (`surface_flux_hourly`, `surface_flux_monthly`):
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

**Solar** (`slv_radiation_hourly`, `slv_radiation_monthly`)
* albedo: surface albedo 
* swgdn: surface incoming shortwave flux  (W m-2)
* swtdn: toa incoming shortwave flux  (W m-2)
* t2m: 2-meter air temperature (K)

**Temperature** (`surface_flux_dailymeans`)
* hournorain: time-during an hour with no precipitation (s)
* tprecmax: Maximum precipitation rate during the period (kg m-2 s-1)
* t2max: max 2-meter air temperature (K)
* t2mmean: mean 2-meter air temperature (K)
* t2mmin: min 2-meter air temperature (K)

**Aerosols** (`surface_aerosol_hourly`)
* bcsmass: black carbon surface mass concentration (kg m-3)
* dusmass25: dust surface mass concentration - pm 2.5 (kg m-3)
* ocsmass: organic carbon surface mass concentration (kg m-3)
* so4smass: SO4 Surface Mass Concentration (kg m-3)
* sssmass25: Sea Salt Surface Mass Concentration - pm 2.5 (kg m-3)

## Outputs

MERRA2 surface flux data currently supports the following outputs:

**Wind** (`surface_flux_hourly`, `surface_flux_monthly`):
* [Wind generation time-series]
* [Wind speed time-series]
* [Wind power density time-series]

**Solar** (`slv_radiation_hourly`, `slv_radiation_monthly`)
* [PV generation time-series]

**Temperature** (`surface_flux_dailymeans`)
* [Celsius Temperature]

**Aerosols** (`surface_aerosol_hourly`)
* [PM2.5 time series]



# Generating Wind Outputs with MERRA2 Data

**geodata** currently supports the following wind outputs using [MERRA2 hourly, single-level surface flux diagnostics](https://disc.gsfc.nasa.gov/datasets/M2T1NXFLX_5.12.4/summary) or also using [MERRA2 monthly mean, single-level surface flux diagnostics](https://disc.gsfc.nasa.gov/datasets/M2TMNXFLX_5.12.4/summary).

* Wind generation time-series (`wind`)
* Wind speed time-series (`windspd`)
* Wind power density time-series (`windpwd`)


## Wind Generation Time-series

Convert wind speeds for turbine to wind energy generation.

```
geodata.convert.wind(
    cutout, 
    turbine, 
    smooth=False, 
    var_height)
```

### Parameters

* `cutout` - **string** -  A cutout created by `geodata.Cutout()`
* `turbine` - **string or dict** - Name of a turbine known by the reatlas client or a turbineconfig dictionary with the keys 'hub_height' for the hub height and 'V', 'POW' defining the power curve.  For a full list of currently supported turbines, see [the list of Turbines here.](https://github.com/east-winds/geodata/tree/master/geodata/resources/windturbine)
* `smooth` - **bool or dict** - If `True`, smooth power curve with a gaussian kernel as determined for the Danish wind fleet to Delta_v = 1.27 and sigma = 2.29. A dict allows to tune these values.

*Note* - 
You can also specify all of the general conversion arguments documented in the `convert_and_aggregate` function (e.g. `var_height='lml'`). 

### Example Code

```
ds_wind = geodata.convert.wind(
                cutout, 
                turbine='Suzlon_S82_1.5_MW', 
                smooth=True, 
                var_height='lml'
            )

ds_wind.to_dataframe(name = 'wind')
```




## Wind Speed Time-series

Extract wind speeds at given height (ms-1)


```
geodata.convert.windspd(
    cutout, 
    **params)
```

### Parameters

* `cutout` - **string** -  A cutout created by `geodata.Cutout()`
* `**params` - Must have 1 of the following:
    - `turbine` - **string or dict** - Name of a turbine known by the reatlas client or a turbineconfig dictionary with the keys 'hub_height' for the hub height and 'V', 'POW' defining the power curve.  For a full list of currently supported turbines, see [the list of Turbines here.](https://github.com/east-winds/geodata/tree/master/geodata/resources/windturbine)
    - `hub-height` - **num** - Extrapolation height (m)
*Note* - 
You can also specify all of the general conversion arguments documented in the `convert_and_aggregate` function (e.g. `var_height='lml'`). 

### Example Code

```
ds_windspd = geodata.convert.windspd(
                cutout, 
                turbine='Vestas_V66_1750kW', 
                var_height='lml'
            )

ds_windspd.to_dataframe(name = 'windspd')
```




## Wind Power Density Time-series

Extract wind power density at given height, according to:
**WPD = 0.5 * Density * Windspd^3**

```
geodata.convert.windwpd(
    cutout, 
    **params)
```

### Parameters

* `cutout` - **string** -  A cutout created by `geodata.Cutout()`
* `**params` - Must have 1 of the following:
    - `turbine` - **string or dict** - Name of a turbine known by the reatlas client or a turbineconfig dictionary with the keys 'hub_height' for the hub height and 'V', 'POW' defining the power curve.  For a full list of currently supported turbines, see [the list of Turbines here.](https://github.com/east-winds/geodata/tree/master/geodata/resources/windturbine)
    - `hub-height` - **num** Extrapolation height (m)
*Note* - 
You can also specify all of the general conversion arguments documented in the `convert_and_aggregate` function (e.g. `var_height='lml'`).

### Example Code

```
ds_windwpd = geodata.convert.windwpd(
                cutout, 
                turbine='Suzlon_S82_1.5_MW', 
                var_height='lml'
            )

ds_windwpd.to_dataframe(name = 'windwpd')
```

# Generating Solar Outputs with MERRA2 Data
**geodata** currently supports the following solar outputs using the following configs:
* `slv_radiation_hourly`: 
    - [MERRA2 hourly, single-level diagnostics](https://disc.gsfc.nasa.gov/datasets/M2SDNXSLV_5.12.4/summary)
    - [MERRA2 hourly, single-level radiation diagnostics](https://disc.gsfc.nasa.gov/datasets/M2T1NXRAD_5.12.4/summary)
* `slv_radiation_monthly`:
    - [MERRA2 monthly mean, single-level diagnostics](https://disc.gsfc.nasa.gov/datasets/M2TMNXSLV_5.12.4/summary)
    - [MERRA2 monthly mean, single-level radiation diagnostics](https://disc.gsfc.nasa.gov/datasets/M2TMNXRAD_5.12.4/summary)

* PV generation time-series (`pv`)

## PV Generation Time-series

Convert downward-shortwave, upward-shortwave radiation flux and
	ambient temperature into a pv generation time-series.

### Parameters

* `cutout` - **string** -  A cutout created by `geodata.Cutout()`
* `panel` - string - Specify a solar panel type on which to base the calculation.  **geodata** contains an internal solar panel dictionary with keys defining several solar panel characteristics used for the time-series calculation.  For a complete list of included panel types, see [the list of panel types here.](https://github.com/east-winds/geodata/tree/master/geodata/resources/solarpanel)
* `orientation` - str, dict or callback - Panel orientation can be chosen from either `latitude_optimal`, a constant orientation such as `{'slope': 0.0,'azimuth': 0.0}`,  or a callback function with the same signature as the callbacks generated by the `geodata.pv.orientation.make_*` functions.
* (optional) clearsky_model - string or None - 	Either the `simple` or the `enhanced` Reindl clearsky model. The default choice of None will choose dependending on data availability, since the `enhanced` model also incorporates ambient air temperature and relative humidity.

### Example Code and Result

```
ds_pv = geodata.convert.pv(
    cutout, 
    panel="KANEKA", 
    orientation = "latitude_optimal"
    )

ds_pv.to_dataframe(name = 'pv')
```

# Generating Aerosol Outputs with MERRA2 Data
**geodata** currently supports the following aerosol output using the following config:
* `surface_aerosol_hourly`: [MERRA2 hourly, single level aerosol diagnostics](https://disc.gsfc.nasa.gov/datasets/M2T1NXAER_5.12.4/summary)

* PM2.5 time series (`pm25`)

## PM2.5 time series

Generate PM2.5 time series according to [1]:
		PM2.5 = [Dust2.5] + [SS2.5] + [BC] + 1.4*[OC] + 1.375*[SO4]

### Parameters

* `cutout` - **string** -  A cutout created by `geodata.Cutout()`

### Example Code and Result

```
ds_pm25 = geodata.convert.pv(
    cutout
    )

ds_pm25.to_dataframe(name = 'pm25')
```

# Example Output

With the full list of supported MERRA2 outputs above, we here provide an example of generating output with [MERRA2 hourly, single-level surface flux diagnostics](https://disc.gsfc.nasa.gov/datasets/M2T1NXFLX_5.12.4/summary).

## Setup

Let's assume we've downloaded MERRA2 data and have created a cutout along the following lines:

```
cutout = geodata.Cutout(name="merra2-europe-sub24-2012-02",
                       module="merra2",
                       weather_data_config="surface_flux_hourly",
                       xs=slice(30, 41.56244222),
                       ys=slice(33.56459975, 35),
                       years=slice(2012, 2012),
                       months=slice(2,2))
cutout.prepare(overwrite=True);
```

We can now use this cutout to generate datasets.

## Creating a wind speed time-series

To create a wind speed time-series, we can use the following code with our MERRA2 cutout:


```
ds = geodata.convert.windspd(cutout, turbine='Suzlon_S82_1.5_MW', var_height='lml')
```

Some information about the parameters:
* `turbine` - string - Specify a turbine type on which to base the calculation.  **geodata** contains an internal turbine dictionary with  eys 'hub_height' for the hub height and 'V' and 'POW' defining the power curve.  For a complete list of included turbine types, see [the list of Turbines here.](https://github.com/east-winds/geodata/tree/master/geodata/resources/windturbine)
* `var_height` - string - suffix for variables containing wind speed and variable height


The convert function returns an xarray dataset, which is an in-memory representation of a NetCDF file:

```
<xarray.DataArray 'wnd79m' (time: 720, y: 3, x: 19)>
array([[[1.6262563 , 1.6435686 , 2.3172426 , ..., 7.8951974 ,
         7.2711506 , 6.5789685 ],
        [1.3850213 , 1.3126061 , 1.9565196 , ..., 5.5188637 ,
         5.551135  , 5.743653  ],
        [1.3309381 , 0.9012143 , 1.3341541 , ..., 3.1683128 ,
         4.024107  , 4.058724  ]]...

Coordinates:
    lon      (x) float64 30.0 30.62 31.25 31.88 32.5 ... 39.38 40.0 40.62 41.25
    lat      (y) float64 34.0 34.5 35.0
  * x        (x) float64 30.0 30.62 31.25 31.88 32.5 ... 39.38 40.0 40.62 41.25
  * y        (y) float64 34.0 34.5 35.0
  * time     (time) datetime64[ns] 2012-11-01T00:30:00 ... 2012-11-30T23:30:00
Attributes:
    long_name:  extrapolated 79.0 m wind speed using log ratio, from variable...
    units:      m s**-1
```

To convert this array to a more conventional dataframe, we can run:

```
df = ds.to_dataframe(name='wind')
```

which converts the xarray dataset into a pandas dataframe. 


The result is a dataset with an observation for each time period (in the MERRA2 case, hourly) and geographic point with the `extrapolated 79.0 m wind speed (units: m s-1)` for each observation (output as indicated by the output xarray dataset metadata).

Finally, we can run something like:

```
df.to_csv('merra_wind_data.csv')
```

to output the data to csv for use in other applications.


