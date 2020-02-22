Readme
-

Builds off library: [atlite](https://github.com/PyPSA/atlite)
See: [README-atlite.rst](README-atlite.rst)


## Installation

### Python


Default:
`python setup.py install`

Local install (use same `prefix` folder as python custom installation):
`python3 setup.py install --prefix=$HOME/opt/python-3.6.9`


**Other required libraries:**
- boto3
- xarray
- netcdf4
- toolz
- pyproj
- scipy
- shapely
- rasterio
- geopandas
- progressbar
- requests
- matplotlib
- dask (`pip3 install dask[complete]`)


### MERRA2 credentials

1. Create a .netrc file in your home directory.
		cd ~ or cd $HOME
		touch .netrc
		echo "machine urs.earthdata.nasa.gov login <uid> password <password>" >> .netrc (where <uid> is your user name and <password> is your Earthdata Login password without the brackets)
		chmod 0600 .netrc (so only you can access it)
2. Create a cookie file. This file will be used to persist sessions across calls to wget or curl.
		cd ~ or cd $HOME
		touch .urs_cookies
		Note: you may need to re-create .urs_cookies in case you have already executed wget without valid authentication.


## Examples

`create_cutout_2.py`
- Create and save a cutout for small sub-region of Europe in one month in `output/*.nc`
- if `module="era5"`: will download raw data to temp folder before extracting cutout

`create_wind.ipynb`
- trims variables of downloaded files if necessary
- creates (or loads existing) cutout files
- example extracts wind power at given height

`get_data.ipynb` / `get_data.py`
- download data routines


