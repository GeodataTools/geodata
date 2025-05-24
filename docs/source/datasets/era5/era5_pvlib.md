# Generating PV Outputs with `pvlib` 

**geodata** also supports solar photovoltaic (PV) generation modeling using the `pvlib` library, integrating solar position, irradiance, and PV system modeling with flexible configuration options.

For more information on `pvlib`, see: [pvlib-python.readthedocs.io](https://pvlib-python.readthedocs.io/en/stable/)

## Supported PV Output


### Solar Photovoltaic Generation Time-series (`pvlib_model`)

Convert ERA5 data into PV generation estimates using `pvlib` functionality.

```python
pvlib_model(
        cutout, 
        system,
        model_chain_config,
        vars = [
            'influx_diffuse', 
            'influx_direct', 
            'dewpoint_temperature',
            'temperature', 
            'wnd100m'
        ]
    ) -> xarray.Dataset
```

#### Parameters

* `cutout` - `xr.Dataset` - An input ERA5 cutout with the following variables:
      - **influx_diffuse** (*float*) - Diffuse horizontal irradiance.  
      - **influx_direct** (*float*) - Direct normal irradiance.  
      - **dewpoint_temperature** (*float*) - Dewpoint temperature in Celsius.  
      - **temperature** (*float*) - Air temperature in Celsius.  
      - **wnd100m** (*float*) - Wind speed at 100m.  
* `system` - A `PVSystem` class defined by `pv_system()`, describing the collection and interactions of PV system components to be used in modeling.
* `model_chain_config` - A `ModelChainConfig` class that defines pvlib ModelChain parameters that can be passed to one or more instances of pvlib_model(). Allows reuse of a common set of ModelChain parameters across multiple PVSystems or even multiple cutouts. 
*Note*: For full information on input parameters for the `PVSystem` and `ModelChainConfig` classes, [see the docstrings here.](https://github.com/GeodataTools/geodata/blob/master/src/geodata/pvlib.py)


#### Output

* Outputs an `xarray.Dataset` containing:
    - **ac** (*float*) - AC photovoltaic output (W).
    - **pv** (*float*) - Photovoltaic capacity.

#### Example Code and Result

```python
n_mods = 50
n_strings = 1
cec_modules = geodata.pvlib.retrieve_sam('CECMod')
module = cec_modules['Kaneka_U_SA105']
inv =  geodata.pvlib.retrieve_sam("CECInverter")['Fronius_USA__CL_33_3_Delta__208V_']

system = geodata.pvlib.pv_system(
    arrays = None,
    surface_tilt=35,
    surface_azimuth=180,
    racking_model = 'open_rack',
    module_parameters=module,
    modules_per_string = n_mods,
    module_type = 'glass_polymer',
    module = 'Kaneka_U_SA105',
    strings_per_inverter = n_strings, 
    inverter_parameters=inv
)

model_config = geodata.pvlib.ModelChainConfig(
    clearsky_model= 'haurwitz',
    transposition_model='perez', 
    solar_position_method= 'nrel_numpy',
    airmass_model= 'kastenyoung1989',
    dc_model='cec',
    ac_model='sandia', 
    aoi_model="physical",
    spectral_model='first_solar',
    dc_ohmic_model='no_loss'
)

model = geodata.pvlib.pvlib_model(
    cutout,
    system,
    model_config
)
model

```
