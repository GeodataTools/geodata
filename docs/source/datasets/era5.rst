ERA5 Specific Instructions
==========================

This page explains how you can set up access to ERA5 data from the `Copernicus Data Store <https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels?tab=overview>`_.

Creating a CDS account
----------------------

To download ERA5 data from the CDS, you'll need to create a free `CDS account here <https://cds.climate.copernicus.eu/>`_.

Download data through CDS API
-----------------------------

Once your account has been created, set up access to the API by following these steps:

1. Log into your CDS account and visit your `profile page <https://cds.climate.copernicus.eu/profile>`_.
2. Install the API key. There will be a section called **Personal Access Token**.
   Copy these two lines into a file called ``.cdsapirc`` in your user root folder.

- **macOS/Linux**: Open a terminal and run:

  .. code-block:: bash

      touch ~/.cdsapirc

  Then add the lines using:

  .. code-block:: bash

      echo [line 1 of the code] >> ~/.cdsapirc
      echo [line 2 of the code] >> ~/.cdsapirc


  - **Windows**: The process is slightly more complicated. Please refer to the in-depth guide at the Copernicus Knowledge Base `here <https://confluence.ecmwf.int/display/CKB/How+to+install+and+use+CDS+API+on+Windows>`_.

3. Install the CDS API client by opening a terminal/shell and running

.. code-block:: bash

     pip install ".[download]"

(Assuming you are in Geodata's *root directory*.)

1. Once you've installed the API key and the API client, confirm access by running an
   example in a Python script or a Jupyter notebook:

.. code-block:: python

     import cdsapi

     c = cdsapi.Client()

     c.retrieve(
          "reanalysis-era5-single-levels",
          {
               "product_type": "reanalysis",
               "format": "netcdf",
               "variable": [
                    "2m_dewpoint_temperature",
                    "2m_temperature",
               ],
               "year": "2011",
               "month": [
                    "01",
               ],
               "day": ["01", "02", "03"],
               "time": [
                    "00:00",
                    "12:00",
               ],
          },
          "download.nc",
     )

The above example downloads 2m temperature and 2m dewpoint temperature with data points
at 00:00 and 12:00 for each day, from January 1-3, 2011, in NetCDF format.

If this works, you have successfully set up access to the ERA5 data through the CDS API.
Please subsequently refer to the `general documentation on datasets <../overview.rst>`_
for more information on how to download ERA5-based datasets using the ``geodata``
package.
