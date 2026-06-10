ERA5 Specific Instructions
==========================

This page covers **CDS account and API credential setup** for ERA5. Once credentials
are in place, use the dataset classes — do not call ``cdsapi`` by hand for routine
downloads.

Recommended download method
-----------------------------

The **recommended way** to fetch ERA5 data in Geodata is:

1. Complete the CDS setup below (one-time).
2. Follow :ref:`downloading-era5-data` in :doc:`overview` — ``load_dataset``,
   instantiate with ``years`` / ``months`` / optional ``bounds``, then ``download()``.

Geodata's ERA5 classes (for example ``ERA5Wind3DHourlyDataset``) create a
``cdsapi.Client`` internally and submit the correct product requests for each
registered ``weather_config``.

Creating a CDS account
----------------------

To download ERA5 data from the CDS, create a free `CDS account here <https://cds.climate.copernicus.eu/>`_.

Configure CDS API credentials
-----------------------------

Once your account exists, install local API access:

1. Log into your CDS account and visit your `profile page <https://cds.climate.copernicus.eu/profile>`_.
2. Under **Personal Access Token**, copy the two lines for your ``.cdsapirc`` file
   (URL and key).

**macOS/Linux** — create ``~/.cdsapirc``:

.. code-block:: bash

   touch ~/.cdsapirc
   # Paste the two lines from your CDS profile into ~/.cdsapirc

**Windows** — see the Copernicus guide on
`installing the CDS API on Windows <https://confluence.ecmwf.int/display/CKB/How+to+install+and+use+CDS+API+on+Windows>`_.

Ensure ``cdsapi`` is available (it is a dependency of Geodata when you install the
package). Then proceed to :ref:`downloading-era5-data` in :doc:`overview`.

Verify CDS API access (optional)
--------------------------------

You can confirm credentials with a minimal ``cdsapi`` script. This is **optional** —
Geodata dataset downloads use the same client and credentials.

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
           "month": ["01"],
           "day": ["01", "02", "03"],
           "time": ["00:00", "12:00"],
       },
       "download.nc",
   )

This example fetches 2 m temperature and dewpoint at 00:00 and 12:00 UTC for
2011-01-01 through 2011-01-03. If it succeeds, your ``.cdsapirc`` is valid.

For production workflows, prefer :ref:`downloading-era5-data` in :doc:`overview` so
Geodata requests the correct ERA5 products, paths, and post-processing for
``wind_3d_hourly``, ``wind_solar_hourly``, and other registered configs.

What's next
-----------

- :ref:`downloading-era5-data` in :doc:`overview` — **recommended** download workflow
- :doc:`../development/offline-era5-fixture-datasets` — offline ``*_test`` configs for CI
- :doc:`../modeling/wind/index` or :doc:`../modeling/pvlib/index` — run models on downloaded data
