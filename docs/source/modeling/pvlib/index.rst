PVLib Modeling
==============

PVLib is a Python library for modeling solar photovoltaic systems. It provides a set of tools for modeling the performance of solar photovoltaic systems

How to use the model
---------------------

The PVLib models are imported from the `pvlib` module.

Step 1: Import the necessary libraries
----------------------------------------

To get started, we need to import the required libraries. We will import
the `pvlib` from the `geodata` library, as well as any other
libraries needed for data handling and visualization.

.. code:: Python

    import xarray as xr

    from geodata.datasets import load_dataset
    from geodata.model.pvlib import Pvlib

Step 2: Load the dataset
------------------------

Next, we need to load the dataset that contains the solar irradiance data.
We will use the `wind_solar_hourly` dataset from the ERA5 dataset.

.. code:: Python

    # Load the dataset
    ds_cls = load_dataset("wind_solar_hourly")
    ds = ds_cls(
        years = slice(2016, 2016),
        months = slice(1, 1)
    )
    if not ds.downloaded:
        ds.download() # Download the data if we don't have it locally
    print(ds.downloaded) # Check if the dataset is downloaded. Should return True.

Step 3: Create the model with specific configs
----------------------------------------------

Next, we need to create the model with specific configs.

.. code:: Python

    model = Pvlib(ds)

Two configurations are required: (1) **PV system setup** — physical array geometry (tilt, azimuth), module and inverter from the SAM database, and racking; (2) **Model config** — algorithms for clearsky irradiance, transposition, solar position, airmass, DC/AC conversion (CEC, Sandia), and losses (AOI, spectral, ohmic).
Following is an example of how to create the model with specific configs.

.. code:: Python

    # create the pv_system
    n_mods = 50
    n_strings = 1
    cec_modules = model.retrieve_sam('CECMod')
    module = cec_modules['Kaneka_U_SA105']
    inv = model.retrieve_sam("CECInverter")['Fronius_USA__CL_33_3_Delta__208V_']
    model.init_pv_system(
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

.. code:: Python

    model.init_model_config(
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

Step 4: Estimate the capacity factor
------------------------------------

Next, we can estimate the AC Power and PV capacity using the model.

.. code:: Python

    cf = model.estimate(
        years = slice(2016, 2016), 
        months = slice(1, 1),
        xs = slice(8, 10), # Optional: longitude subset
        ys = slice(48, 46), # Optional: latitude subset (see below)
    )
    print(cf)

The output will be an xarray Dataset containing the estimated AC power (``ac``) and
capacity factor (``pv``) for the specified region and time period.

Estimate options
----------------

All models inherit a common pattern for **time** and **space** subsetting via
``estimate()``. The PVLib model adds one extra output option.

Temporal subsetting
~~~~~~~~~~~~~~~~~~~

Pass ``years`` and ``months`` as ``slice`` objects to limit the period processed.
Omit either argument to use the prepared model's full range (subject to what was
available when ``prepare()`` ran).

Spatial subsetting (``xs``, ``ys``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Pass ``xs`` and ``ys`` as ``slice(start, stop)`` to restrict longitude (``x``) and
latitude (``y``). Omit either argument to keep the full horizontal extent of the
prepared dataset.

Geodata **normalizes slice bounds** before calling xarray's ``.sel()``. You can pass
bounds in either order (for example ``ys=slice(48, 46)`` for a band in central
Europe) and still get a non-empty selection. This matters for ERA5-style grids where
latitude is often stored in **descending** order: a naive ``slice(46, 48)`` would
return no points without normalization.

.. code:: Python

    # Equivalent selections on a descending-latitude grid:
    cf_a = model.estimate(years=slice(2016, 2016), months=slice(1, 1), ys=slice(48, 46))
    cf_b = model.estimate(years=slice(2016, 2016), months=slice(1, 1), ys=slice(46, 48))

Compact output (``compact_output``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

By default, ``estimate()`` returns a compact dataset with only two data variables:

- ``ac`` — AC power (W)
- ``pv`` — capacity factor (AC output normalized by module nameplate)

Set ``compact_output=False`` to retain **all intermediate weather and ModelChain
columns** per grid cell (irradiance components, temperature, wind, and other inputs
used along the chain). Use this for debugging or when you need columns beyond
``ac`` and ``pv``; the result is larger and slower to write.

.. code:: Python

    # Default: only ac and pv
    cf = model.estimate(
        years=slice(2016, 2016),
        months=slice(1, 1),
        compact_output=True,
    )
    list(cf.data_vars)  # ['ac', 'pv']

    # Full per-coordinate table (debugging / downstream analysis)
    full = model.estimate(
        years=slice(2016, 2016),
        months=slice(1, 1),
        compact_output=False,
    )
    list(full.data_vars)  # ac, pv, plus weather and intermediate columns

.. toctree::
   :maxdepth: 1
   :caption: Tutorials on Specific Models

