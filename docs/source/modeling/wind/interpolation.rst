Tutorial: Estimate Wind Speed with Interpolation
================================================

In this tutorial, we will learn how to estimate wind speed using the interpolation model from the geodata library.

.. warning::
   Performing wind speed estimation using interpolation requires a dataset with known
   wind speed values at **multiple** locations.

   Currently, only the :code:`weather_data_config` :code:`wind_3d_hourly` from the ERA5 dataset
   contains the necessary wind speed data for interpolation.

   Therefore, all of the information below only applies with :code:`wind_3d_hourly` or cutouts
   derived from it. Using any other dataset will lead to a :code:`ValueError`.

Step 1: Import the necessary libraries
----------------------------------------

To get started, we need to import the required libraries. We will import the `WindInterpolationModel` from the `geodata` library, as well as any other libraries needed for data handling and visualization.

.. code:: Python

    import geodata
    import xarray as xr

    from geodata.model.wind import WindInterpolationModel

Step 2: Load the dataset
------------------------

Next, we need to load the dataset that contains the wind speed data. We will use the `wind_3d_hourly` dataset from the ERA5 dataset.

.. code:: Python

    # Load the dataset
    ds = geodata.Dataset(
        module="era5",
        weather_data_config="wind_3d_hourly",
        years=slice(2006, 2006),
        months=slice(1, 1),
        bounds=[-10, 35, 10, 45] # Optional: specify the bounding box
    )

    if not ds.prepared:
        ds.get_data()  # Download the data if we don't have it locally

Step 3: Prepare the data for interpolation
--------------------------------------------
Before we can perform the actual estimation process, we need to compute the parameters
required for the interpolation model.

.. code:: Python

    # Create a model based on the above dataset. The model will be associated with
    # the dataset forever. If you wish to use a different dataset, you will need to
    # create a new model.

    model = WindInterpolationModel(ds)
    model.prepare()

If you have already prepared a cutout with :code:`wind_3d_hourly`, you can also pass
that into the model as well. The model treats dataset and cutouts indifferently.
Simply replace :code:`ds` with your cutout variable.

.. note::
   The `prepare` method computes the necessary parameters for the interpolation model
   based on the loaded dataset. Everything will be saved under the :code:`models`
   directory under the path :code:`GEODATA_ROOT`.

.. note::
    It is not necessary to call the `prepare` method every time you want to perform
    interpolation. You only need to call it once after loading the dataset. From that
    point on, you can load and use the model directly without re-preparing it.

Step 4: Perform the interpolation
----------------------------------

Now that we have prepared the model, we can perform the interpolation to estimate wind
speed at the desired locations. Suppose we want to estimate the wind speed at a height
of 60 above ground during Jaunary of 2006 for the entire region covered by the original
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

.. note::
    As the underlying ERA5 dataset already contained wind speed at certain heights, the
    model also has a feature to return the original wind speed values from the dataset
    if desired. To do this, simply set the `use_real_data` parameter to `True` in the
    `estimate` method. You do not need worry about whether the height you queried is
    available in the dataset; the model will handle that for you. If the height is not
    available, it will perform interpolation instead.
