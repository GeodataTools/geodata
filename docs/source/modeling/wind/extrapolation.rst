Tutorial: Estimate Wind Speed with Extrapolation
================================================

In this tutorial, we will learn how to estimate wind speed using the extrapolation model
 from the geodata library.

.. warning::
   Performing wind speed estimation using extrapolation requires a dataset with known
   wind speed values at **multiple** locations.

   Currently, only the :code:`weather_data_config` :code:`slv_flux_hourly` from the MERRA2 dataset
   contains the necessary wind speed data for extrapolation.

   Therefore, all of the information below only applies with :code:`slv_flux_hourly` or cutouts
   derived from it. Using any other dataset will lead to a :code:`ValueError`.

Step 1: Import the necessary libraries
----------------------------------------

To get started, we need to import the required libraries. We will import the `WindExtrapolationModel` from the `geodata` library, as well as any other libraries needed for data handling and visualization.

.. code:: Python

    import xarray as xr

    from geodata.datasets import load_dataset
    from geodata.model.wind import WindExtrapolationModel


Step 2: Load the dataset
------------------------

Next, we need to load the dataset that contains the wind speed data. We will use the `slv_flux_hourly` dataset from the ERA5 dataset.

.. code:: Python

    # Load the dataset
    ds_cls = load_dataset("slv_flux_hourly")
    ds = ds_cls(
        years=slice(2006, 2006),
        months=slice(1, 1),
        bounds=[-10, 35, 10, 45]  # Optional: specify the bounding box
    )

    if not ds.downloaded:
        ds.download()  # Download the data if we don't have it locally

    print(ds.downloaded)  # Check if the dataset is downloaded. Should return True.


Step 3: Compute extrapolation parameters
--------------------------------------------
The extrapolation is separated into two steps, first estimating extrapolation parameters
using linear regression, and second extrapolating to desired heights.
First, we compute the extrapolation parameters.
For more information on the model, see the section below: `How the Extrapolation Model Works`_.

.. code:: Python

    # Create a model based on the above dataset. The model will be associated with
    # the dataset forever. If you wish to use a different dataset, you will need to
    # create a new model.

    model = WindExtrapolationModel(ds)
    model.prepare()

If you have already prepared a cutout with the config :code:`slv_flux_hourly`, you
can also pass
that into the model as well. The model treats dataset and cutouts indifferently.
Simply replace :code:`ds` with your cutout variable.

.. note::
   The `prepare` method computes the necessary parameters for the extrapolation model
   based on the loaded dataset. Everything will be saved under the :code:`models`
   directory under the path :code:`GEODATA_ROOT`.

.. note::
    It is not necessary to call the `prepare` method every time you want to perform
    extrapolation. You only need to call it once after loading the dataset. From that
    point on, you can load and use the model directly without re-preparing it.

Step 4: Estimate using the extrapolation model
----------------------------------------------

Now that we have prepared the model, we can perform the extrapolation to estimate wind
speed at the desired locations. Suppose we want to estimate the wind speed at a height
of 60 above ground during January of 2006 for the entire region covered by the original
dataset, we can do this as follows:

.. code:: Python

    estimated_wind_speed = model.estimate(
        height=60,
        years=slice(2006, 2006),
        months=slice(1, 1),
    )

This will return an xarray DataArray containing the estimated wind speed values. Note
that you can also select a subset area by passing in :code:`xs=slice(start, end)`
and/or :code:`ys=slice(start, end)` parameters to the `estimate` method.


Step 5: Estimate Wind Turbine Capacity Factor (CF) using the interpolation model
--------------------------------------------------------------------------------

Geodata also supports a limited set of wind turbine models to estimate the capacity
factor (CF) of a wind turbine directly. To get a list of available wind turbine models,
you can use the `get_available_windturbines` function:

.. code:: Python

    from geodata.resource import get_available_windturbines

    turbines = get_available_windturbines()
    print(turbines)  # List of available wind turbine configurations


To estimate the capacity factor of a wind turbine, you can use the `estimate` method
and passign in the `turbine` parameter with the name of the wind turbine model.

.. code:: Python

    # Estimate the capacity factor for a specific wind turbine model
    estimated_cf: xr.Dataset = model.estimate(
        turbine="Vestas_V112_3MW",  # Example wind turbine model
        years=slice(2006, 2006),
        months=slice(1, 1),
    )

    print(estimated_cf)  # Display the estimated capacity factor


The output will be an xarray Dataset containing the estimated capacity factor values
for the specified wind turbine model over the given time period and region.


How the Extrapolation Model Works
---------------------------------

The model calculates hub height wind speed from MERRA2, extrapolating the variables in
MERRA's tavg1_2d_slv_Nx data collection, which is a set of the time-averaged
single-layer diagnostics.

Specifically, the variables we use for extrapolation are: 2-m wind (U2M, V2M, in m/s),
10-m wind (U10M, V10M), 50-m wind (U50M, V50M), and the zero-plane displacement
height (DISPH, in meters). Additionally, we also use the wind speed at MERRA2's lowest
model level (ULML, VLML, in m/s), the height of the lowest model level
(HLML, in meters), may vary depending on the location. We can obtain the wind speed at
any given location and height by computing the norm of the vector sum of the U and V
components.


The hub height wind speed can be calculated as

.. math::
   \nu  = \alpha \ln\left(\frac{H - d}{z}\right)

.. math::
   z = e^{-\beta/\alpha}

where :math:`\nu` is the hub height wind speed, :math:`\alpha` is the best-fit slope
from a linear regression of wind speeds on vertical heights, :math:`\ln` is the natural logarithm, :math:`H` is the hub height,
:math:`d` is the zero-plane displacement height, and :math:`\beta` is the intercept
from the linear regression fit.

Here we estimate :math:`\alpha` and :math:`\beta` fitting a simple linear regression model to the heights and wind speeds in the data.
