# get_data.py
# 	Builds wind profile from cutout (eg create_cutout_2.py)

import logging
logging.basicConfig(level=logging.INFO)

import atlite

# import matplotlib.pyplot as plt
# outdir = 'output/'

DS = atlite.Dataset(	module="merra2",
						xs=slice(30, 41.56244222),
						ys=slice(40, 33.56459975),
						years=slice(2011, 2011),
						months=slice(1,1))
print(DS)
if DS.prepared == False:
	DS.get_data()

# DS.set_saved_files()

DS.trim_variables()
