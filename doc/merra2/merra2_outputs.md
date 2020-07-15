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
    panel="KANENA", 
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