Back to the [Table of Contents](https://github.com/GeodataTools/geodata/blob/master/doc/general/tableofcontents.md).

# Modelling Wind Speed Using Geodata With MERRA2 Datasets

Geodata offers a intuitive method for wind speed estimation using MERRA2 datasets. Below is an walkthrough on how such task is done.

## Dataset Preparation

Currently, Geodata only supports dataset objects. The support for dataset-based cutout objects will be added in the future.

```python
# Optional, execute these lines if you want more verbose logs.
import logging

logging.basicConfig(level=logging.INFO)
```

We will use the `slv_flux_hourly` weather config. The time period will be between January 2010 and March 2010.
```python
import geodata

dataset =  geodata.Dataset(
    module="merra2",
    weather_data_config="slv_flux_hourly",
    years=slice(2010, 2010),
    months=slice(1, 3)
)

if not dataset.prepared:
    dataset.get_data()
```

## Compute Model Parameter

For this example, we will use extrapolation model. Simply put, given a series of heights and associated wind speeds, wind speed at other heights could be estimated with: 
$$\nu = \alpha \ln \left( \frac{H - d}{\exp(-\beta / \alpha)} \right)$$

where $\nu$ is the hub height wind speed, $H$ is the hub height, $d$ is the zero-plane displacement height. When computing the parameters, the model would use the measured wind speed at three levels (2m, 10m, 50m, and surface level height) to compute the best-fit coefficient $\alpha$ and intercept $\beta$. 

```python
model = geodata.model.WindModel(dataset)

if not model.prepared:
    model.prepare()
```

A progress bar will be display to indicate the computation progress. After prepration, we are ready to estimate the wind speed at any height, in any geographic region.

Once the model is prepared, all computed parameters will be locally saved in the directory specified by `GEODATA_ROOT`. 

## Estimate Wind Speed

With a prepared model we could estimate wind speed with the following method:

```python
speed = model.estimate(
    height=42, 
    years=slice(2010, 2010), 
    months=slice(1, 1)
    xs=slice(67, 125),
    ys=slice(25, 50)
)
```

This will used the stored parameters to estimate windspeed throughout the continental United States during January 2010. In this example, similar to other geodata modules, we specified the longitude and latitude ranges through `xs` and `ys` argument.

### Optional Parameters

It is worth noting that, with the exception of `height` and `years`, all other arguments are optional. For instance, if a user does not specify `months`, then all months (January through December) will be estimated. The same applies with `xs` and `ys` arguments.

### Result

Currently, the model module would output a xarray `DataArray` object, containing wind speeds indexed through spatial and temporal coordiantes. For more information on visualizing wind speeds with this `DataArray`, please refer to the [xarray documentation](https://docs.xarray.dev/en/stable/) for more guidance.

In the future, we tentatively plan to use a `ModelResult` module as the represnetation for model estimation results. This should allow better interoperability between the model module and other modules of the library.
