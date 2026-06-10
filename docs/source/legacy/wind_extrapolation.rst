Wind extrapolation (legacy)
===========================

.. note::

   **Legacy / untested in CI.** ``WindExtrapolationModel`` only supports the
   ``slv_flux_hourly`` weather config (MERRA-2 via ``load_dataset``). It is **not**
   part of the current ERA5 workflow documented on the homepage. For ERA5 wind, use
   :doc:`/modeling/wind/interpolation` instead.

For the recommended modern path, see the :doc:`documentation homepage </index>`.

Tutorial: Estimate Wind Speed with Extrapolation
------------------------------------------------

In this tutorial, we will learn how to estimate wind speed using the extrapolation model
from the geodata library.

.. warning::

   Extrapolation requires a dataset with wind speed at **multiple** heights. In Geodata,
   only ``slv_flux_hourly`` (MERRA-2) is registered for ``WindExtrapolationModel``.
   Using any other ``weather_config`` raises ``ValueError``.

Step 1: Import the necessary libraries
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: Python

    import xarray as xr

    from geodata.datasets import load_dataset
    from geodata.model.wind import WindExtrapolationModel


Step 2: Load the dataset
~~~~~~~~~~~~~~~~~~~~~~~~

Use the MERRA-2 ``slv_flux_hourly`` config (not ERA5):

.. code:: Python

    ds_cls = load_dataset("slv_flux_hourly")
    ds = ds_cls(
        years=slice(2006, 2006),
        months=slice(1, 1),
        bounds=[-10, 35, 10, 45],
    )

    if not ds.downloaded:
        ds.download()

    print(ds.downloaded)


Step 3: Compute extrapolation parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: Python

    model = WindExtrapolationModel(ds)
    model.prepare()

Prepared coefficients are stored under ``GEODATA_ROOT/models/`` (see
:doc:`/modeling/wind/index` — **Preparing the model** for ``prepare`` / ``prepared`` /
``force``).

Step 4: Estimate using the extrapolation model
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: Python

    estimated_wind_speed = model.estimate(
        height=60,
        years=slice(2006, 2006),
        months=slice(1, 1),
    )

Step 5: Estimate wind turbine capacity factor
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

See :doc:`/modeling/wind/interpolation` Step 5 — **Understanding the output** for what
``cf`` means. Example:

.. code:: Python

    estimated_cf = model.estimate(
        turbine="Vestas_V112_3MW",
        years=slice(2006, 2006),
        months=slice(1, 1),
    )


How the Extrapolation Model Works
---------------------------------

The model calculates hub height wind speed from MERRA2 surface and low-level winds,
extrapolating variables in MERRA's tavg1_2d_slv_Nx collection (2 m, 10 m, 50 m winds,
displacement height, lowest model level winds, etc.).

The hub height wind speed uses a log-profile fit (see the original tutorial in the
repository history for the full equations).
