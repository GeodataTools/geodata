Back to the [Table of Contents](https://github.com/east-winds/geodata/blob/master/doc/general/tableofcontents.md).

Wind parameterizations
-

_To read/write latex formula, install GitHub add-on: [TeXify](https://github.com/apps/texify)._

## Extrapolation

Extrapolating wind speeds from measured (or modeled) heights to other heights must assume a model of the atmosphere. This is particularly important for studying wind power generation as turbine hub heights can be 50-120m above the ground. Standard approaches include the log law and power law.

Geodata allows users to customize the wind extrapolation function through the `extrap_fn` parameter in the call to `wind.extrapolate_wind_speed()`. By default, the log law ratio (`extrap_fn = log_ratio`) is used.


## Stability correction

Simple extrapolation routines assume neutral stability of the atmosphere, potentially ignoring highly stable or unstable atmospheric conditions that could alter the wind speed profile at certain hours of the day.


**`L` Obukhov length**

Most stability correction functions rely on the Obukhov length, which can be calculated from available parameters [[1]](#1):

<img src="/doc/parameterizations/tex/85b8a1a1ba5b0d126cb83dc1fc238bbc.svg?invert_in_darkmode&sanitize=true" align=middle width=116.15572649999999pt height=50.332382100000004pt/>
<img src="/doc/parameterizations/tex/42c4a538b0898977509d72e9843ffaad.svg?invert_in_darkmode&sanitize=true" align=middle width=185.66718059999997pt height=50.332382100000004pt/>
<img src="/doc/parameterizations/tex/46f5f4c88affca88279076c9fd0b8502.svg?invert_in_darkmode&sanitize=true" align=middle width=103.41612599999999pt height=28.091038800000003pt/>

where <img src="/doc/parameterizations/tex/142a0e55e7be4656275f11bc455886d4.svg?invert_in_darkmode&sanitize=true" align=middle width=16.145467799999988pt height=14.15524440000002pt/> is the friction velocity,  <img src="/doc/parameterizations/tex/acbc3da08b00b7616f45697dcc62d44d.svg?invert_in_darkmode&sanitize=true" align=middle width=15.527067599999992pt height=28.091038800000003pt/> is the virtual temperature, <img src="/doc/parameterizations/tex/bd971e27e71d9bf4792553ddfe23d537.svg?invert_in_darkmode&sanitize=true" align=middle width=15.62926694999999pt height=14.15524440000002pt/> is the density, <img src="/doc/parameterizations/tex/7db7bb668a45f8c25263aa8d42c92f0e.svg?invert_in_darkmode&sanitize=true" align=middle width=18.52532714999999pt height=22.465723500000017pt/> is the specific heat, <img src="/doc/parameterizations/tex/c39017bed62ee8b2a521779ae976ebe0.svg?invert_in_darkmode&sanitize=true" align=middle width=52.39338884999999pt height=21.18721440000001pt/>, and <img src="/doc/parameterizations/tex/547e731741446703f67244e558fc0508.svg?invert_in_darkmode&sanitize=true" align=middle width=27.205208249999988pt height=22.465723500000017pt/> is the virtual heat flux in terms of sensible and latent heat fluxes (negative if directed upwards).


**<img src="/doc/parameterizations/tex/191de8bb417c3c87a59bd9e2f4222ed0.svg?invert_in_darkmode&sanitize=true" align=middle width=22.372663499999987pt height=22.831056599999986pt/> Stability correction**

Extrapolated wind speeds are "corrected" via a stability correction function <img src="/doc/parameterizations/tex/191de8bb417c3c87a59bd9e2f4222ed0.svg?invert_in_darkmode&sanitize=true" align=middle width=22.372663499999987pt height=22.831056599999986pt/> according to:

<img src="/doc/parameterizations/tex/b215101f20292f1e8b0eb0ee28cf1d61.svg?invert_in_darkmode&sanitize=true" align=middle width=171.50142899999997pt height=24.7161288pt/>

where <img src="/doc/parameterizations/tex/5421d2d8ce8d3d6a7eb5c5e6c23527a5.svg?invert_in_darkmode&sanitize=true" align=middle width=63.754082699999984pt height=24.65753399999998pt/> takes different forms in the literature [[2]](#2).

Geodata has the following stability correction functions installed (in terms of the parameter <img src="/doc/parameterizations/tex/084c4d336ce7b83477dc39e2638fa8c5.svg?invert_in_darkmode&sanitize=true" align=middle width=27.774070499999986pt height=24.65753399999998pt/>):
- `psi_linear`: linear for positive parameters (This function does not perform well for high <img src="/doc/parameterizations/tex/084c4d336ce7b83477dc39e2638fa8c5.svg?invert_in_darkmode&sanitize=true" align=middle width=27.774070499999986pt height=24.65753399999998pt/>)
- `psi_linearexp`: piecewise linear-exponential
- `psi_linearexpconst`: piecewise linear-exponential with maximum constant correction, illustrated here for various constants:
<img src="wind_stability_corrections.png" width="400">


## References

<a id="1">[1]</a> Rose, S., & Apt, J. (2016). Quantifying sources of uncertainty in reanalysis derived wind speed. Renewable Energy, 94, 157–165. https://doi.org/10.1016/j.renene.2016.03.028

<a id="2">[2]</a> Sharan, M., & Aditi. (2009). Performance of various similarity functions for nondimensional wind and temperature profiles in the surface layer in stable conditions. Atmospheric Research, 94(2), 246–253. https://doi.org/10.1016/j.atmosres.2009.05.014
