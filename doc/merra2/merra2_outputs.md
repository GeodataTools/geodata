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
* `turbine` - **string or dict** - Name of a turbine known by the reatlas client or a turbineconfig dictionary with the keys 'hub_height' for the hub height and 'V', 'POW' defining the power curve.  For a full list of currently supported turbines, see [LINK]
* `smooth` - **bool or dict** - If `True`, smooth power curve with a gaussian kernel as determined for the Danish wind fleet to Delta_v = 1.27 and sigma = 2.29. A dict allows to tune these values.

*Note* - 
You can also specify all of the general conversion arguments documented in the `convert_and_aggregate` function (e.g. `var_height='lml'`). See [LINK]

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
    - `turbine` - **string or dict** - Name of a turbine known by the reatlas client or a turbineconfig dictionary with the keys 'hub_height' for the hub height and 'V', 'POW' defining the power curve.  For a full list of currently supported turbines, see [LINK]
    - `hub-height` - **num** - Extrapolation height (m)
*Note* - 
You can also specify all of the general conversion arguments documented in the `convert_and_aggregate` function (e.g. `var_height='lml'`). See [LINK]

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
    - `turbine` - **string or dict** - Name of a turbine known by the reatlas client or a turbineconfig dictionary with the keys 'hub_height' for the hub height and 'V', 'POW' defining the power curve.  For a full list of currently supported turbines, see [LINK]
    - `hub-height` - **num** Extrapolation height (m)
*Note* - 
You can also specify all of the general conversion arguments documented in the `convert_and_aggregate` function (e.g. `var_height='lml'`). See [LINK]

### Example Code

```
ds_windwpd = geodata.convert.windwpd(
                cutout, 
                turbine='Suzlon_S82_1.5_MW', 
                var_height='lml'
            )

ds_windwpd.to_dataframe(name = 'windwpd')
```