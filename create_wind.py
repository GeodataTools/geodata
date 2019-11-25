# create_wind.py
# 	Builds wind profile from cutout (eg create_cutout_2.py)

import logging
logging.basicConfig(level=logging.INFO)

import atlite
import matplotlib.pyplot as plt

outdir = 'output/'

## Load cutout
# if already created, should just return reference to that folder
cutout = atlite.Cutout(name="europe-sub-2011-01",
                       module="era5",
                       xs=slice(30, 41.56244222),
                       ys=slice(40, 33.56459975),
                       years=slice(2011, 2011),
                       months=slice(1,1,1) )
# cutout.prepare()

## Wind profiles
# call: wind(cutout, turbine, smooth=False, **params)
#	.. cutout.convert_and_aggregate(convert_func=convert_wind, turbine=turbine,**params)
# Returns xr.DataArray
ds = atlite.convert.wind(cutout, turbine='Suzlon_S82_1.5_MW', smooth=True)

# Pandas and save to csv
df = ds.to_dataframe(name='power')
df.to_csv(outdir + 'era5_wind_test.csv')

# Plot
df.plot.scatter(x='x', y='y', c='power') #, cmap='jet',vmin=0, vmax=250)
plt.colorbar()
plt.show()
