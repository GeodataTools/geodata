# get_data.py
# 	Download dataset for given module, years, months

import logging
logging.basicConfig(level=logging.INFO)

import atlite

DS = atlite.Dataset(	module="merra2",
						years=slice(2011, 2011),
						months=slice(1,1))
print(DS)
if DS.prepared == False:
	DS.get_data()

# set all saved files to modify via trim (as opposed to recently downloaded)
# DS.set_saved_files()

# eliminate unneeded variables (to save filespace)
DS.trim_variables()
