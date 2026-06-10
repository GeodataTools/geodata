# ERA5 model outputs

After you download ERA5 data (see :ref:`downloading-era5-data` in the
[Dataset module overview](../datasets/overview.rst)), **models** turn raw files into
analysis-ready time series. Masking applies afterward on model results (see
[Mask tutorials](../mask/xarray_mask_tutorial.ipynb)).

This page is a short catalog of common **model** outputs on the current tested path. It
does not describe legacy ``Cutout`` / ``geodata.convert`` products (see
[Legacy MERRA2 outputs](../legacy/merra2/merra2_outputs.md)).

## Wind generation time-series

Hub-height **capacity factor** (``cf``) from ERA5 3D wind:

| Step | Component |
|------|-----------|
| Dataset | ``wind_3d_hourly`` (or ``wind_3d_hourly_test`` for offline fixtures) |
| Model | ``WindInterpolationModel`` — [Wind modeling](wind/index.rst), [interpolation tutorial](wind/interpolation.rst) |

```python
from geodata.datasets import load_dataset
from geodata.model.wind import WindInterpolationModel

ds = load_dataset("wind_3d_hourly")(years=slice(2016, 2016), months=slice(1, 1))
if not ds.downloaded:
    ds.download()

model = WindInterpolationModel(ds)
model.prepare()
cf = model.estimate(turbine="Vestas_V112_3MW", years=slice(2016, 2016), months=slice(1, 1))
```

## Wind speed time-series

Same dataset and model; pass ``height=<meters>`` instead of ``turbine=``:

```python
wind_speed = model.estimate(height=100.0, years=slice(2016, 2016), months=slice(1, 1))
```

## Solar photovoltaic generation time-series

Hourly **AC power** (``ac``) and **capacity factor** (``pv``):

| Step | Component |
|------|-----------|
| Dataset | ``wind_solar_hourly`` (or ``wind_solar_hourly_test`` for offline fixtures) |
| Model | ``Pvlib`` — [PVLib modeling](pvlib/index.rst) |

After ``init_pv_system()`` and ``init_model_config()``, call ``estimate()`` (see the
pvlib docs for ``compact_output`` and spatial subsetting).

## See also

- [ERA5 CDS setup](../datasets/era5.rst)
- [Offline ERA5 fixtures](../development/offline-era5-fixture-datasets.md)
- [Supported input/output formats](../quick_start/input_output.md)
