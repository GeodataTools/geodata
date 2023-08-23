# ERA5 Dataset

This page explains how you can setup access to ERA5 data from the [Copernicus Data Store](https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels?tab=overview).

## Creating a CDS account

To download the ERA5 data from the CDS, you'll need to create a free [CDS account here](https://cds.climate.copernicus.eu/user/register).

## Download data through CDS API
Once your account has been created, you'll need to set up access to the API by doing the following:

1.  Log into your CDS account and visit [this page](https://cds.climate.copernicus.eu/api-how-to).  
2.  Install the API key.  The uppermost code block on the right side will display the 
    api url and your user key.  You will then need to copy these two lines this into 
    a file called _.cdsapirc_ that you will create in your user root folder.
     - *macOS/Linux*: For macOS or Linux, you can create _.cdsapirc_ by opening a terminal window and running `touch ~/.cdsapirc` to create the file and then `echo [line 1 of the code] >> ~/.cdsapirc` followed by `echo [line 2 of the code] >> ~/.cdsapirc`.
     - *Windows*: For Windows, the process is slightly more complicated. Please refer to an in-depth guide at the Copernicus Knowledge Base [here](https://confluence.ecmwf.int/display/CKB/How+to+install+and+use+CDS+API+on+Windows).
3.  Install the CDS API client by opening a terminal/shell and running `pip install ".[download]"`, assuming you are in Geodata's *root directory*.
4.  Once you've installed the API key and the API client, you can confirm access by running an example like those from the [CDS API how-to](https://cds.climate.copernicus.eu/api-how-to) in a Python script or a Jupyter notebook.

```python
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
```

The above example downloads 2m temperature and 2m dewpoint temperature with data points at 00:00 and 12:00 for each day, from January 1-3, 2011, in NetCDF format.

## More Resources on ERA5

We've also provide a few tutorials on how to interact with the ERA5 datasets using Geodata. They can be found in the links below.

```{toctree}
:maxdepth: 2
:glob:
*
```
