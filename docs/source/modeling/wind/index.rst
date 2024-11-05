Wind Modeling
=============

Starting from geodata v0.2.0, geodata's capability to model and estimate wind speed have
been from the cutout module to a separate wind module. This module has the capability to
estimate wind speed with two modes: interpolation and extrapolation.

How to use the models
---------------------

Unlike other modules of geodata, the model module does not get imported automatically.
In other words, one cannot use the models with an import statement like this:

.. code:: Python

    import geodata
    model = geodata.model.wind.WindInterpolationModel()

Instead, the user must import the models explicitly:

.. code:: Python

    from geodata.model.wind import WindInterpolationModel
    model = WindInterpolationModel()

The reason for this is to keep the main geodata namespace clean,
since we might add many more models in the future.

More details on how to use the models can be found in each model's respective tutorial
as well as in the API reference.

.. toctree::
   :maxdepth: 1
   :caption: Tutorials on Specific Models

   interpolation
   extrapolation
