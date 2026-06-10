Tutorial: Estimate Wind Speed with Interpolation
================================================

In this tutorial, we will learn how to estimate wind speed using the interpolation model
 from the geodata library.

.. warning::
   Performing wind speed estimation using interpolation requires a dataset with known
   wind speed values at **multiple** locations.

   Currently, only the :code:`weather_data_config` :code:`wind_3d_hourly` from the ERA5 dataset
   contains the necessary wind speed data for interpolation.

   Therefore, all of the information below only applies with :code:`wind_3d_hourly` or cutouts
   derived from it. Using any other dataset will lead to a :code:`ValueError`.

Step 1: Import the necessary libraries
----------------------------------------

To get started, we need to import the required libraries. We will import
the `WindInterpolationModel` from the `geodata` library, as well as any other
libraries needed for data handling and visualization.


.. code:: Python

    import xarray as xr

    from geodata.datasets import load_dataset
    from geodata.model.wind import WindInterpolationModel


Step 2: Load the dataset
------------------------

Next, we need to load the dataset that contains the wind speed data.
We will use the `wind_3d_hourly` dataset from the ERA5 dataset.


.. code:: Python

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


Step 3: Compute interpolation parameters
--------------------------------------------
The interpolation is separated into two steps to separate the computationally-intensive
step (estimating interpolation parameters) from the computationally-easy step
(interpolating at desired heights).

.. note::
    Wind-speed estimation using **cubic spline interpolation** fits a smooth, piecewise cubic polynomial
    to known wind-speed data across spatial or temporal dimensions. Given a set of data points
    :math:`(x_0, y_0), (x_1, y_1), \dots, (x_n, y_n)`, where :math:`x_i` are known positions (e.g., time or altitude)
    and :math:`y_i` are corresponding wind speeds, the cubic spline for each interval
    :math:`[x_i, x_{i+1}]` is defined as:

    .. math::

        S_i(x) = a_i + b_i(x - x_i) + c_i(x - x_i)^2 + d_i(x - x_i)^3

    The coefficients :math:`a_i, b_i, c_i, d_i` are determined by solving a system of equations subject to:

    1. **Interpolation condition**:
    :math:`S_i(x_i) = y_i`, and :math:`S_i(x_{i+1}) = y_{i+1}`
    2. **Continuity of first derivative**:
    :math:`S_i'(x_{i+1}) = S_{i+1}'(x_{i+1})`
    3. **Continuity of second derivative**:
    :math:`S_i''(x_{i+1}) = S_{i+1}''(x_{i+1})`
    4. **Boundary conditions**, typically *natural*:
    :math:`S_0''(x_0) = 0`, and :math:`S_{n-1}''(x_n) = 0`

    The resulting spline provides a smooth and continuous estimate of wind speed, allowing accurate
    interpolation between measured data points. With `wind_3d_hourly`, we have wind speeds at seven different heights AGL.
    This enables us to estimate wind speeds at any height AGL within that range (from approximately 10m AGL to 170m AGL).

First, we compute the interpolation parameters.

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

Step 4: Estimate wind speeds using the interpolation model
----------------------------------------------------------

Now that we have prepared the model, we can perform the interpolation to estimate wind
speed at the desired locations. Suppose we want to estimate the wind speed at a height
of 60 m above ground during January of 2006 for the entire region covered
by the original dataset, we can do this as follows:

.. code:: Python

    estimated_wind_speed: xr.Dataset = model.estimate(
        height=60.0,
        years=slice(2006, 2006),
        months=slice(1, 1),
    )


This will return an xarray Dataset containing the estimated wind speed values. You can
also restrict the horizontal domain with ``xs`` and ``ys`` (see
:doc:`/modeling/wind/index` — **Estimate options** for slice-order behavior on
ERA5 grids).

.. code:: Python

    estimated_wind_speed = model.estimate(
        height=60.0,
        years=slice(2006, 2006),
        months=slice(1, 1),
        xs=slice(8, 10),
        ys=slice(48, 46),
    )

Step 5: Estimate Wind Turbine Capacity Factor (CF) using the interpolation model
--------------------------------------------------------------------------------

Geodata also supports a limited set of wind turbine models to estimate the capacity
factor (CF) of a wind turbine directly. To get a list of available wind turbine models,
you can use the ``get_available_windturbines`` function:

.. code:: Python

    from geodata.resource import get_available_windturbines

    turbines = get_available_windturbines()
    print(turbines)  # List of available wind turbine configurations


Pass the YAML **stem** (filename without ``.yaml``) as ``turbine`` — for example
``Vestas_V112_3MW`` for ``src/geodata/resources/windturbine/Vestas_V112_3MW.yaml``.

.. code:: Python

    estimated_cf = model.estimate(
        turbine="Vestas_V112_3MW",
        years=slice(2006, 2006),
        months=slice(1, 1),
    )

    print(estimated_cf)

Understanding the output
~~~~~~~~~~~~~~~~~~~~~~

``estimate(turbine=...)`` returns an ``xarray.DataArray`` named ``cf`` with dimensions
``(time, x, y)`` when those coordinates are present.

Geodata computes CF in three steps:

1. **Hub-height wind speed** — interpolate to the turbine's ``HUB_HEIGHT`` from the
   YAML (same vertical spline as Step 4, but at the turbine height rather than a
   height you pass manually).
2. **Power from the power curve** — map wind speed to power (MW) by interpolating the
   tabulated ``V`` / ``POW`` pairs in the turbine YAML.
3. **Normalize** — ``cf = power / P``, where ``P`` is the rated power (maximum value
   in ``POW``).

So ``cf`` is a **dimensionless capacity factor** in ``[0, 1]`` (values can exceed 1
briefly if the curve extrapolates above rated power). Values outside the tabulated
wind-speed range use SciPy's ``interp1d`` extrapolation — treat edge cases with care
in sensitivity analysis.

For implementation details, see ``WindBaseModel._estimate_power`` in the
:ref:`API reference <modindex>`.
