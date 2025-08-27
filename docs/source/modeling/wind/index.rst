Wind Modeling
=============

Starting from geodata v0.2.0, geodata's capability to model and estimate wind speed have
been from the cutout module to a separate wind module. This module has the capability to
estimate wind speed with two modes: interpolation and extrapolation.

How to use the models
---------------------

Unlike some other modules of geodata, the model module does not get imported
automatically. In other words, one cannot use the models with an import statement
like this:

.. code:: Python

    import geodata
    model = geodata.model.wind.WindInterpolationModel()

Instead, the user must import the models explicitly:

.. code:: Python

    from geodata.model.wind import WindInterpolationModel
    model = WindInterpolationModel()

The reason for this is to keep the main geodata namespace clean,
since we might add many more models in the future.

In general, models are created based off of a dataset object, we must be downloaded
first. We'll use the `wind_3d_hourly` dataset from the ERA5 dataset as an example:

.. code:: Python

    from geodata.datasets import load_dataset
    from geodata.model.wind import WindInterpolationModel

    # Load the dataset
    ds_cls = load_dataset("wind_3d_hourly")
    ds = ds_cls(
        years=slice(2006, 2006),
        months=slice(1, 1),
        bounds=[-10, 35, 10, 45]  # Optional: specify the bounding box
    )

    if not ds.downloaded:
        ds.download()  # Download the data if we don't have it locally

    print(ds.downloaded)  # Check if the dataset is downloaded. Should return True.


Once we have the dataset, we can create a model based on it.
The model will be associated with the dataset forever, so if you wish to use a
different dataset, you will need to create a new model.

.. code:: Python

    # Create a model based on the above dataset
    model = WindInterpolationModel(ds)
    model.prepare()  # Prepare the model

    print(model.prepared) # Check if the model is prepared. Should return True.


Once the model is prepared, we can use it to estimate wind speed at desired heights.

.. code:: Python

    # Estimate wind speed at a specific height (e.g., 100 meters)
    wind_speed: xr.Dataset = model.estimate(height=100.0)

    # The wind_speed variable is an xarray Dataset containing the estimated wind speed at the specified height.
    # Over the region covered by the original dataset.


The above demonstrates the typical workflow. More model-specific details can be found
in each model's respective tutorial as well as in the API reference.

.. toctree::
   :maxdepth: 1
   :caption: Tutorials on Specific Models

   interpolation
   extrapolation
