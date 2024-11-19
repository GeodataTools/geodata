Wind Models
=============

Estimating wind speed and wind generation from simulated reanalysis or weather models entails a number of assumptions, due to the general sparsity of measurements and complex meteorology affecting wind speed evolution with height. 
Default approaches use parameterizations such as the log law to scale wind speed to desired turbine hub heights according to general relationships. 
The `model` module extracts multiple wind speeds from the underlying datasets which are used to estimate wind speed in one of two ways: interpolation and extrapolation.

How to use the wind models
---------------------

Unlike other modules of geodata, the `model` module does not get imported automatically, in order to keep the namespace clean. 

Instead, the user must import `model` explicitly:

.. code:: Python

    from geodata.model.wind import WindInterpolationModel
    model = WindInterpolationModel()

For more details on how to use the models, see the tutorials below
as well as the API reference.

.. toctree::
   :maxdepth: 1
   :caption: Tutorials on Specific Models

   interpolation
   extrapolation
