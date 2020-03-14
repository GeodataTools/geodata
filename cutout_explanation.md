



### ERA5

#### Fixed errors

There were two errors, which I fixed in this commit:
- `cutout.py`, 176-177
- `era5.py`, 277-279

Both MERRA2 and ERA5 now work for me. I demonstrate them both in `cutout_test.ipynb`.

More generally, I thought I would explain a little the issue you noticed with `orography`, which relates to the structure of the cutout process, which is sort of convoluted, but is designed to speed up some operations. I include below some ideas for improving ERA5's representation, in particular.

#### 1. Metadata -- Download dimensions of files via metadata

Idea: download a small file that has all of the metadata of the eventual (larger) files in order to determine the necessary dimensions / coordinates (in `xarray` parlance)

This is located in `preparation.py`:
```
def cutout_get_meta:
	...
	meta_kwds = cutout.meta_data_config.copy()
	meta_kwds.update(dataset_params)

	# Get metadata (eg prepare_meta_merra2)
	prepare_func = meta_kwds.pop('prepare_func')
```
Here, `prepare_func` comes from `cutout.meta_data_config`
(in case of ERA5: `meta_data_config = dict(prepare_func=prepare_meta_era5)`)

Then, if you look at the `prepare_meta_era5` definition:
```
def prepare_meta_era5(xs, ys, year, month, module):
	# Load/download ERA5 height data as metadata
	...
	with _get_data(variable='orography',
```

#### 2. Download data

This occurs in `era5.py`, in the function `prepare_month_era5`, which calls:
```
	with _get_data(area=area, year=year, month=month,
				   variable=[
					   '100m_u_component_of_wind',
				...
```


```
def _get_data(target=None, product='reanalysis-era5-single-levels', chunks=None, **updates):
	"""
	Check local folders for ERA5
	If necessary: download ERA5 data from the Climate Data Store (CDS)
	"""

	## Check if local copy
	# filename = first characters of each variable
	f = "".join([v[0] for v in updates['variable']])
	if target is None:
		target = era5_dir

	fn=os.path.join(era5_dir, '{year}/{month:0>2}/{f}.nc'.format(year=updates['year'], month=updates['month'],f=f))

	if os.path.isfile(fn):
	#	Local file exists
		with xr.open_dataset(fn, chunks=chunks) as ds:
		...
	else:
	#	Download new file
		...
```
Here, because multiple variables are downloaded to a single file, I rename the file to just include the first letter of each variable name (e.g., `112rsssstt.nc` for the wind and solar variables).

The calculation then proceeds using these files.

#### Some options to improve

1. Create an ERA5 dataset download step (like with MERRA2). This would involve recoding the `_get_data()` function, such that one could call something like:
```
DS = atlite.Dataset(module="era5",
					years=slice(2011, 2011),
					months=slice(1,1))
DS.get_data()
```

2. Check for coverage in the downloaded files, and if they do not contain the desired cutouts, return a warning. Provide tips for downloading the larger file (via the current method using the cutout and API)

3. Something else I haven't thought of...


### MERRA2

#### Data

Here are the global files I have downloaded (2 days of the month). Note these are hourly:
```
gps-mrdavidson:01 michd$ pwd
/Users/michd/Documents/GEODATA/data/merra2/2011/01
gps-mrdavidson:01 michd$ ls -l
total 458792
-rw-------  1 michd  staff  114823776 Nov 26 18:43 MERRA2_400.tavg1_2d_flx_Nx.20110101.nc4
-rw-------  1 michd  staff  115431457 Nov  5 11:53 MERRA2_400.tavg1_2d_flx_Nx.20110102.nc4
```

These are different than the one you downloaded, which is monthly means (see https://cmr.earthdata.nasa.gov/search/concepts/C1276812865-GES_DISC.html):
`MERRA2_400.tavgM_2d_adg_Nx.201101.nc4`

You can use the `Dataset` class to download the appropriate hourly files (e.g., see: `get_data.ipynb`). The hourly files are indicated in `merra2.py`:
```
weather_data_config = {
	...
	'surface_flux': dict(...
						template=os.path.join(merra2_dir, '{year}/{month:0>2}/MERRA2_*.tavg1_2d_flx_Nx.*.nc4'),
						url = 'https://goldsmr4.gesdisc.eosdis.nasa.gov/data/MERRA2/M2T1NXFLX.5.12.4/{year}/{month:0>2}/MERRA2_{spinup}.tavg1_2d_flx_Nx.{year}{month:0>2}{day:0>2}.nc4',
						fn = os.path.join(merra2_dir, '{year}/{month:0>2}/MERRA2_{spinup}.tavg1_2d_flx_Nx.{year}{month:0>2}{day:0>2}.nc4'),
						...
						)
```

To create a monthly average, I think it makes sense to add a new element in `weather_data_config` (e.g., `'surface_flux_monthly'`). For now, you might try to get the hourly data working first.

#### Cutout

This script now works for me:
```
import logging
logging.basicConfig(level=logging.INFO)
import atlite

## Load cutout
# if already created, should just return reference to that folder
cutout = atlite.Cutout(name="merra2-europe-sub2-2011-01",
                       module="merra2",
                       xs=slice(30, 41.56244222),
                        ys=slice(33.56459975, 35),
                       years=slice(2011, 2011),
                       months=slice(1,1,1) )
cutout.prepare();
```

Generating:
```
gps-mrdavidson:merra2-europe-sub2-2011-01 michd$ ls -l
total 448
-rw-r--r--@ 1 michd  staff  195706 Mar 13 20:32 201101.nc
-rw-r--r--@ 1 michd  staff   29227 Mar 13 20:32 meta.nc
```


