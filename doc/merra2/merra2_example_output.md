# MERRA2 Example Output

A short guide to generating output with [MERRA2 hourly, single-level surface flux diagnostics](https://disc.gsfc.nasa.gov/datasets/M2T1NXFLX_5.12.4/summary).

For a full list of supported MERRA2 outputs, see: INSERT LINK HERE

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
* `turbine` - string - Specify a turbine type on which to base the calculation.  **geodata** contains an internal turbine dictionary with  eys 'hub_height' for the hub height and 'V' and 'POW' defining the power curve.  For a complete list of included turbine types, see [the list of Turbines here.](INSERT LINK HERE)
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

which converts the xarray dataset into a pandas dataframe.  Let's examine the output:

[image link]


The result is a dataset with an observation for each time period (in the MERRA2 case, hourly) and geographic point with the `extrapolated 79.0 m wind speed (units: m s-1)` for each observation (output as indicated by the output xarray dataset metadata).

Finally, we can run something like:

```
df.to_csv('merra_wind_data.csv')
```

to output the data to csv for use in other applications.
