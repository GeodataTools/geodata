# Mask on cutout workflow

A short guide to use the product of mask module for masking cutout objects.

For a complete jupyter notebook example, see: [mask_on_cutout_test.ipynb](https://github.com/east-winds/geodata/tree/mask/tests/mask_on_cutout_test.ipynb)

## 1. Introduction

Geodata is able to process geospatial data to extract cutouts over specified geographies. Built off the [rasterio library](https://rasterio.readthedocs.io/en/latest/quickstart.html), the mask module imports rasters and shapefiles, merges and flattens multiple layers together, and extracts subsetted cutout data from merged masks and shapefiles.

The sample workflow that shows how to create a mask object and save it to the disk can be found in `mask_creation_workflow.md` and [mask_test.ipynb](https://github.com/east-winds/geodata/tree/mask/tests/mask_test.ipynb).

## 2. Setup

To start, import the required dependencies:

```
import geodata
import xarray as xr
```

To launch a logger for detailed debugging, run:

```
import logging
logging.basicConfig(level=logging.INFO)
```

We will use a cutout object created through:

```
geodata.Dataset(module="merra2", years=slice(2011, 2011),month=slice(1,12),weather_data_config = "slv_radiation_monthly")
```

```
geodata.Cutout(name = "china-2011-slv-test",module = "merra2",weather_data_config = "slv_radiation_monthly", xs = slice(73, 136),ys = slice(18, 54),years = slice(2011, 2011), months = slice(1,12))
```

## 3. Loading mask

The code block below shows how the china_bin mask is created. Once the code below runs successfully, th user would not need to re-run the code everytime.

```
prov_path = shpreader.natural_earth(resolution='10m', category='cultural', name = 'admin_1_states_provinces')
china_all_shapes = geodata.mask.get_shape(prov_path, key = 'name_en', return_dict = True,
                         condition_key = 'admin', condition_value = 'China')

#Remove the islands province, it does not have a key in the shape dictionary. Therefore we would remove it by its key 'None'.
china_all_shapes.pop(None) 

china_bin = geodata.Mask("china_bin")
china_bin.add_layer('FINAL_GRID_5BINS.tif', layer_name = 'bins')

#Extract shape on bins layer
china_bin.extract_shapes(china_all_shapes, layer = 'bins')
china_bin.save_mask()
```

If the user have already created such mask object on disk, it is sufficient to run the following cell to retrieve the object:

```
china_bin = geodata.mask.load_mask("china_bin")
```

## 4. Adding mask variables to the cutout object

### 4.1 Adding masking variables

`add_mask` method will add attribute `merged_mask` and `shape_mask` from the mask object to the cutout object. However, the `merged_mask` or `shape_mask` in the cutout object will be stored in the format of xarray.DataArray, and their dimensions will be coarsened to the same dimension with the cutout metadata.

The mask `china_bin` has no `merged_mask` value, but the `add_mask` method will look for both `merged_mask` and `shape_mask` attribute saved for the loaded mask, unless the user set the parameter `merged_mask = False`, or `shape_mask = False`

```
cutout.add_mask("china_bin")
```

### 4.2 Adding area variable

The user is also able to add grid area for each grid cell in the cutout metadata. Because the grid cell with the same latitude difference will have different area due to the cylindrical map projection, this method makes sure that the user can capture the variation of grid cell area in the dataset.

```
cutout.add_grid_area()
```

### 4.3 Creating PV data through cutout conversion

The code block below will use the `geodata.convert.pv` method to generate `ds_cutout`, an xarray Dataset that contains the pv variable for the cutout.

```
ds_solar = geodata.convert.pv(cutout, panel = "KANENA", orientation = "latitude_optimal")
ds_solar = ds_solar.reset_coords(['lon', 'lat'], drop = True)
ds_solar = ds_solar.rename({'x': 'lon', 'y': 'lat'})
ds_cutout = ds_solar.to_dataset(name = 'solar')
ds_cutout = ds_cutout.coarsen(time = 12, boundary = 'exact').mean()
ds_cutout = ds_cutout.transpose("time", "lat", "lon")
```

### 4.4 Combining variables

The `mask` method will mask converted dataSet variable, such as `ds_cutout` created above, with previously added masks and area variable, and return a dictionary of xarray Dataset. Each key in the dictionary is one unique mask from either the `merged_mask` or `shape_mask` variable from the cutout object, and each value is an xarray dataset containing the dataSet variable (`ds_cutout`) with the mask and area values.

The program will automatically search for `merged_mask` and `shape_mask`, unless the user specify `merged_mask = False` or `shape_mask = False`, the masks in `shape_mask` will have the same key as it has in the `shape_mask` attribute, and the mask for `merged_mask` will have the key name "merged_mask".

```
>>> combine = cutout.mask(dataset = ds_cutout)
>>> combine.keys()
dict_keys(['Anhui', 'Beijing', 'Chongqing', 'Fujian', 'Gansu', 'Guangdong', 'Guangxi Zhuang Autonomous Region', 'Guizhou', 'Hainan', 'Hebei', 'Heilongjiang', 'Henan', 'Hubei', 'Hunan', 'Inner Mongolia', 'Jiangsu', 'Jiangxi', 'Jilin', 'Liaoning', 'Ningxia Hui Autonomous Region', 'Qinghai', 'Shaanxi', 'Shandong', 'Shanghai', 'Shanxi', 'Sichuan', 'Tianjin', 'Tibet', 'Xinjiang', 'Yunnan', 'Zhejiang'])
```

Below we can show the three variables in the 'Anhui' dataset in the `combined` dictionary we generate.

```
combine['Anhui']['solar'].plot()
```

![anhui_slope](https://github.com/east-winds/geodata/blob/mask/images/mask_on_cutout_test/anhui_solar.png)

```
combine['Anhui']['mask'].plot()
```

![anhui_mask](https://github.com/east-winds/geodata/blob/mask/images/mask_on_cutout_test/anhui_mask.png)

```
combine['Anhui']['area'].plot()
```

![anhui_area](https://github.com/east-winds/geodata/blob/mask/images/mask_on_cutout_test/anhui_area.png)
