Readme
-


## Installation

### Python

**Library: `atlite`** (See `README-atlite.rst`)

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
- ipython[complete]
- jupyter (ipython notebook no longer supported)

### MERRA

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


## To Do

1. Download data routines
	- √ era5: `_get_data` in `era5`
  - (ncep is six-hourly: SKIP)
2. √ Add MERRA, MERRA2 dataset download functions
	- download authorization (https://disc.gsfc.nasa.gov/data-access#mac_linux_wget)
	- download and extract only variables needed (`trim_variables`)
4. √ Cutout extension to MERRA
 - ? Rethink `cutout` concept (avoid data duplication)
5. √ Additional wind speed extrapolation functions (see `matlab/merra`)
6. √ in `re-validation/` repo → Extract cells of given sites and years (`matlab/merra/GD_sites.txt`), format ~ `GD_merra`
7. √ Add stability correction functions
	- `wind_corrections.ipynb`: stability corrections `stabilityCorrection_linear`, `stabilityCorrection_linearexp`
	- Obukhov lengths



## Issues

### `L` Obukhov lengths

**Old**
$L= - \dfrac{u_{* }^{3}\overline{\theta_{v}}\rho_{a}c_{p,d}}{\kappa gH_{f}}$

for `a=z/L, z=80m`:
![](output/a_example.png)

**New**

Rose & Apt SI has more complete equation for `L`:
$L= - \dfrac{u_{* }^{3}\overline{\theta_{v}}\rho_{a}c_{p,d}}{\kappa gH_{v0}}$
$H_{v0}=H_{f}+0.61C_{p}\dfrac{\overline{\theta}}{L_{e}}H_{L}$
$\overline{\theta}=\overline{T}\left(p_{0}/p\right)^{\kappa_{p}}$

involves several additional parameters
- Implemented, but not much difference: virtual heat flux (incl latent heat flux, potential temp)
- IGNORE FOR NOW: virtual temp (includes humidity)
	> The difference between the actual and the virtual temperature is small for cold air masses and low specific humidity, but can be several degrees for warm and very humid air masses. (Emeis, p. 19)
- IGNORE FOR NOW: in terms of $C_p = C_{pd} (1+0.84*q)$ , not $C_{pd}$

### $\psi_m$ Stability correction

`wind_corrections.ipynb`
- various stability correction fns
- `_linear`, `_linearexp`, `_linearexpconst`
- default: `linearexpconst`
  - linear, exponential and plateau
  - const = 7 = x/L (above, which $\psi$ is constant)

linearexp:
![](output/wsp_flux_linearexp_example.png)
(var = log_law, var2 = log_law_linearexp)


### xarray, dask dependencies (#1)

```
Traceback (most recent call last):
...
AttributeError: module 'dask' has no attribute 'base'
```

**Solution**
Running `pytest` in `xarray` folder (`/usr/local/lib/python3.7/site-packages/xarray`)

Need module `toolz`

### xarray, dask dependencies (#2)

```
AttributeError: module 'dask' has no attribute 'multiprocessing'
```

**Solution**
`pip3 install dask[complete]`
(https://docs.dask.org/en/latest/install.html#pip)

### ProgressBar

```
Traceback (most recent call last):
...
  File "/Users/michd/Dropbox (MIT)/git/geodata/atlite/convert.py", line 149, in convert_and_aggregate
    maybe_progressbar = make_optional_progressbar(show_progress, prefix, len(yearmonths))
  File "/Users/michd/Dropbox (MIT)/git/geodata/atlite/utils.py", line 38, in make_optional_progressbar
    maybe_progressbar = pgb.ProgressBar(prefix=prefix, widgets=widgets, max_value=max_value)
TypeError: __init__() got an unexpected keyword argument 'prefix'
```

**(Non-)Solution**
Ignore progressbar, set to not display, in `convert.py`:
`
maybe_progressbar = make_optional_progressbar(False, prefix, len(yearmonths))
`

### Long process time during cutout
at step: `Merging variables into monthly compound files`
?
