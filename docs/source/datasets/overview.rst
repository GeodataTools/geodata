==========================
Dataset Module Overview
==========================

The ``geodata.datasets`` module provides tools and classes for accessing, managing, and
processing geospatial datasets. It offers a unified interface for loading various
data formats, handling metadata, and performing common geospatial operations.

Key Features
------------

- Supports the download and management of datasets from various sources, such as
  `ERA5 <https://www.ecmwf.int/en/forecasts/dataset/ecmwf-reanalysis-v5>`_ and
  `MERRA2 <https://gmao.gsfc.nasa.gov/reanalysis/MERRA-2/>`_.

- Provides a consistent API for accessing geospatial data, regardless of the underlying
  data source.

Typical Usage
-------------

In the following example, we will demonstrate how to download a dataset containing wind
and solar data from ECMWF's ERA5 dataset.

.. code-block:: python

    from geodata.datasets import load_dataset

    dataset_cls = load_dataset("wind_solar_hourly")

    years = slice(2010, 2020)
    months = slice(1, 13)
    dataset = dataset_cls(years=years, months=months)

Here, we first create a dataset class using the `load_dataset` function, specifying the
name of the dataset we want to load. We then instantiate the dataset class with the
desired time range (years and months). Then, we can create a dataset instance with
that class, which will handle the downloading and processing of the data.

Dataset Classes
-----------------
The `geodata.datasets` module includes several dataset classes, each tailored for
specific datasets. These classes encapsulate the logic for downloading, processing, and
accessing the data. Some of the available dataset classes
(listed by `weather_data_config`) include:

- `wind_solar_hourly`: A dataset containing hourly wind and solar data from ECMWF's
  ERA5. It is important to note that the wind data are only recorded at
  10 and 100 meters above ground level. Hence, this dataset is also referred to as
  2D wind and solar dataset.

- `wind_3d_hourly`: A dataset containing hourly wind data from ECMWF's ERA5 at
  multiple vertical levels, providing a more comprehensive view of the wind profile.
  It can be used for wind speed estimation using and interpolation model built into
  the geodata library.

You can use the `list_datasets` function to see all available datasets in the
`geodata.datasets` module. This function returns a list of dataset names that can be
loaded using the `load_dataset` function. For example:

.. code-block:: python

    from geodata.datasets import list_datasets

    available_datasets = list_datasets()
    print(available_datasets)  # Outputs a list of available dataset names.

Check Preparedness of Datasets
------------------------------------------------
To check if a dataset is prepared and ready for use, you can use the `downloaded`
property of the dataset instance. This property returns a boolean indicating whether the
dataset is fully prepared. If the dataset is not prepared, you can call the `prepare`
method to download and process the data. For example:

.. code-block:: python

    print(dataset.downloaded)  # Check if the dataset is downloaded. Outputs False here.

    if not dataset.downloaded:
        dataset.download()

    print(dataset.downloaded)  # Outputs True after downloading.

Dataset's Interoperability with Cutout
------------------------------------------------

At the moment, the dataset classes are not interoperable with the `Cutout` class.
In the future, we plan to consolidate the functionalities of the `Cutout` class into the
dataset classes and the modeling module (see :doc:`here<../modeling/wind/index>`).

For now, after downloading a dataset, a good point to move forward would be to use the
:doc:`modeling module <../modeling/index>` to create a model that can do certain types
of modeling with the dataset, such as wind speed estimation or solar PV generation.
