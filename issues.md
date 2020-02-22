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
<img src="/tex/75173e569baa9a8789b2c565023aa533.svg?invert_in_darkmode&sanitize=true" align=middle width=122.26786109999998pt height=50.332382100000004pt/>

for `a=z/L, z=80m`:
![](output/a_example.png)

**New**

Rose & Apt SI has more complete equation for `L`:
<img src="/tex/ece36c7399f34962a3b28aad0c68b7f0.svg?invert_in_darkmode&sanitize=true" align=middle width=122.26786109999998pt height=50.332382100000004pt/>
<img src="/tex/42c4a538b0898977509d72e9843ffaad.svg?invert_in_darkmode&sanitize=true" align=middle width=185.66718059999997pt height=50.332382100000004pt/>
<img src="/tex/46f5f4c88affca88279076c9fd0b8502.svg?invert_in_darkmode&sanitize=true" align=middle width=103.41612599999999pt height=28.091038800000003pt/>

involves several additional parameters
- Implemented, but not much difference: virtual heat flux (incl latent heat flux, potential temp)
- IGNORE FOR NOW: virtual temp (includes humidity)
	> The difference between the actual and the virtual temperature is small for cold air masses and low specific humidity, but can be several degrees for warm and very humid air masses. (Emeis, p. 19)
- IGNORE FOR NOW: in terms of <img src="/tex/932ed822fe69d84187c423ef190b3e5b.svg?invert_in_darkmode&sanitize=true" align=middle width=161.22789869999997pt height=24.65753399999998pt/> , not <img src="/tex/d60774afb1f5f2bb975cc84b078668a5.svg?invert_in_darkmode&sanitize=true" align=middle width=25.36840514999999pt height=22.465723500000017pt/>

### <img src="/tex/191de8bb417c3c87a59bd9e2f4222ed0.svg?invert_in_darkmode&sanitize=true" align=middle width=22.372663499999987pt height=22.831056599999986pt/> Stability correction

`wind_corrections.ipynb`
- various stability correction fns
- `_linear`, `_linearexp`, `_linearexpconst`
- default: `linearexpconst`
  - linear, exponential and plateau
  - const = 7 = x/L (above, which <img src="/tex/7e3c241c2dec821bd6c6fbd314fe4762.svg?invert_in_darkmode&sanitize=true" align=middle width=11.29760774999999pt height=22.831056599999986pt/> is constant)

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