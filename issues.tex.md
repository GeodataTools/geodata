Issue tracking
-

_To read/write latex formula, install GitHub add-on: [TeXify](https://github.com/apps/texify)._

### √ Download data routines

- √ era5: `_get_data` in `era5`
- (ncep is six-hourly: SKIP)

### √ Add MERRA, MERRA2 dataset download functions
- download authorization (https://disc.gsfc.nasa.gov/data-access#mac_linux_wget)
- download and extract only variables needed (`trim_variables`)

### √ Cutout extension to MERRA
- ? Rethink `cutout` concept (avoid data duplication)
- √ Additional wind speed extrapolation functions (see `matlab/merra`)

### √ in `re-validation/` repo
- Extract cells of given sites and years (`matlab/merra/GD_sites.txt`), format ~ `GD_merra`

### √ Add stability correction functions
- `wind_corrections.ipynb`: stability corrections `stabilityCorrection_linear`, `stabilityCorrection_linearexp`
- Obukhov lengths


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