# Incorporating Mask into Cutout Workflow

## 1. Introduction

Geodata is able to process geospatial data to extract cutouts over specified geographies. Built off the [rasterio library](https://rasterio.readthedocs.io/en/latest/quickstart.html), the **mask** module imports rasters and shapefiles, merges and flattens multiple layers together, and extracts subsetted cutout data from merged masks and shapefiles.

After we create a mask, we can incorporate the suitability mask object/file into the Cutout. The cutouts are subsets of data based on specific time and geographic ranges. For more information on the creation of cutout, refer to these tutorials: [Creating Cutouts with MERRA2 Data](https://github.com/east-winds/geodata/blob/master/doc/merra2/merra2_createcutout.md), [Downloading and Creating Cutouts with ERA5 Data](https://github.com/east-winds/geodata/blob/master/doc/era5/era5_download.md)

## 2. Setup

To start, import the geodata package and required libraries.

```python
import geodata
import xarray as xr
import matplotlib.pyplot as plt
```

To launch a logger for detailed debugging, run:

```python
import logging
logging.basicConfig(level=logging.INFO)
```

### 2.1 Download data

We will use a Cutout object created from a downloaded dataset. **If you have already created a cutout, load it here and skip to step 3.**

We first download the dataset through `geodata.Dataset()`. In `get_data()`, if we specify `testing = true`, the program downloads only first file in download list (e.g., first day of month)

```python
DS_hourly_test_chn = geodata.Dataset(module = "merra2", 
                                 years = slice(2011, 2011),
                                 months = slice(1,1),
                                 weather_data_config = "slv_radiation_hourly")
```

```python
if DS_hourly_test_chn.prepared == False: 
    DS_hourly_test_chn.get_data(testing = True)
  
# DS_hourly_test_chn.trim_variables()
```

Extract the cutout from the trimmed dataset.

```python
cutout = geodata.Cutout(name = "china-2011-slv-hourly-test",
                        module = "merra2",
                        weather_data_config = "slv_radiation_hourly",
                        xs = slice(73, 136), 
                        ys = slice(18, 54), 
                        years = slice(2011, 2011), 
                        months = slice(1,1))
cutout.prepare()
```

True
## 3. Load Mask

In this tutorial, we use the `china` mask, created in this documentation: [mask_creation_workflow](https://github.com/east-winds/geodata/blob/master/doc/mask/mask_creation_workflow.md)

```python
# View the contents of the china mask
geodata.mask.load_mask("china")
```
Mask china: 

4 layers: ['elevation_filtered', 'modis_filtered', 'protected', 'slope_filtered'].

Merged_mask merged/flattened. 

3 shape_mask: ['Jiangsu', 'Shanghai', 'Zhejiang']. 

Mask has been saved. 


## 4. Adding Mask variables to the Cutout object

### 4.1 Adding masking variables

The `add_mask` method will add attribute `merged_mask` and `shape_mask` from the Mask object to the Cutout object. Once the mask is added to the Cutout object, the `merged_mask` or `shape_mask` from the Mask object will be stored in the format of xarray.DataArray in the Cutout object, and their dimensions will be coarsened to the same dimension with the Cutout metadata.

The `add_mask` method will look for both `merged_mask` and `shape_mask` attribute saved for the loaded mask, unless the user set the parameter `merged_mask = False`, or `shape_mask = False`.

```python
cutout.add_mask("china")
```
Plot the merged mask, coarsened to cutout resolution

```python
cutout.merged_mask.plot()
```
<matplotlib.collections.QuadMesh at 0x1ee15b33a90>
![png](https://github.com/east-winds/geodata/blob/master/images/mask_on_cutout_workflow/output_21_1.png)

### 4.2 Adding area variable

To calculate and add the variation of grid cell areas by latitude to the cutout, use the `add_grid_area` method. Keeping track of the area for each grid cell is necessary for analyses such as calculating the weighted sum of the grid cells based on their area.

```python
cutout.add_grid_area()
```
### 4.3 Creating PV data through cutout conversion

The code block below will use the `geodata.convert.pv` method to generate `ds_cutout`, an xarray Dataset that contains the pv variable for the cutout. We will also reset the variables names such as `x` and `y` to `lon` and `lat`.

```python
ds_solar = geodata.convert.pv(cutout, panel = "KANEKA", orientation = "latitude_optimal")
ds_solar = ds_solar.reset_coords(['lon', 'lat'], drop = True)
ds_solar = ds_solar.rename({'x': 'lon', 'y': 'lat'})
len(ds_solar.time)
```
24


Next, we transform the xarray DataArray into a xarray DataSet (which can contain multiple DataArray). We also need to remove the time dimension.

Here, we also calculate daily means via `ds_cutout.coarsen(time = 24, boundary = 'exact').mean()`, which aggregates the values over its 24 timestamps.

```python
ds_cutout = ds_solar.to_dataset(name = 'solar')
ds_cutout_mean = ds_cutout.coarsen(time = 24, boundary = 'exact').mean()
```
### 4.4 Combining PV data with Mask

The `mask` method for the Cutout will mask converted xarray.Dataset variable, such as `ds_cutout` and `ds_cutout_mean` created above, by combining it with merged_mask or shape_mask in the Cutout object. It will return a dictionary of xarray Dataset. Each key in the dictionary is one unique mask from either the merged_mask or shape_mask variable from the Cutout object, and each value is an xarray dataset containing the dataSet variable (`ds_cutout` or `ds_cutout_mean`) with the mask and area values.

The program will automatically search for `merged_mask` and `shape_mask` to combine with the xarray.Dataset, unless the user specify `merged_mask = False` or `shape_mask = False`. The masks in `shape_mask` will have the same key as it has in the `shape_mask` attribute, and the mask for `merged_mask` will have the same key name `merged_mask`, as `merged_mask` is unique to each mask.

**Daily averaged PV values**

```python
combine_mean = cutout.mask(dataset = ds_cutout_mean)
combine_mean.keys()
```
dict_keys(['merged_mask', 'Jiangsu', 'Shanghai', 'Zhejiang'])
From the output variable `combine_mean`, check out the combined xarray.Dataset for the Jiangsu province, and plot each of its xarray.DataArray.

```python
combine_mean['Jiangsu']
```
![png](https://github.com/east-winds/geodata/blob/master/images/mask_on_cutout_workflow/output_37_1.png)

Visualize the averaged PV value for each grid cell in the Cutout. Note that the data is the aggregated value for the date.

```python
combine_mean['Jiangsu']['solar'].plot()
```
<matplotlib.collections.QuadMesh at 0x1ee15dd8a90>

![png](https://github.com/east-winds/geodata/blob/master/images/mask_on_cutout_workflow/xarray_output1.png)

Visualize the masking value for each grid cell in the Cutout.

```python
combine_mean['Jiangsu']['mask'].plot()
```
<matplotlib.collections.QuadMesh at 0x1ee15e72520>
![png](https://github.com/east-winds/geodata/blob/master/images/mask_on_cutout_workflow/output_39_1.png)

**Area- and mask-weighted hourly PV values**

We use the raw hourly output generated by cutout to create time-series PV plots weighted by the mask and area. Note that we transposed ds_cutout so that time is set as the first dimension, which ease the following calculation since we want to aggregate the array spatially from each grid cell.

```python
combine = cutout.mask(dataset = ds_cutout.transpose("time", "lat", "lon"))
combine.keys()
```
dict_keys(['merged_mask', 'Jiangsu', 'Shanghai', 'Zhejiang'])


Calculate the aggregated mean solar PV for each provinces, at each time point. We will apply this equation below to calculate the area-weighted average. We save the result into a dictionary `PV_dict`, where its keys are the provinces, and the corresponding values are the PV series.

![png](https://github.com/east-winds/geodata/blob/master/images/mask_on_cutout_workflow/equation.png)

```python
PV_dict = {}

for prov_name in list(combine.keys())[1:]:

    PV_dict[prov_name] = (((combine[prov_name]['solar'] * 
                          combine[prov_name]['mask'] * 
                          combine[prov_name]['area']
                           ).sum(axis = 1).sum(axis = 1)) / 
                         (combine[prov_name]['mask'] * combine[prov_name]['area']).sum())
```
The aggregated PV time-series for Zhejiang province.

```python
PV_dict['Zhejiang']
```
![png](https://github.com/east-winds/geodata/blob/master/images/mask_on_cutout_workflow/xarray_output2.png)

Finally, for each province, plot the solar series weighted by mask * area.

```python
for prov_name, series in PV_dict.items():
  
    plt.plot(series, label = prov_name)
  
    plt.title(f"Solar series weighted by area for Chinese provinces.")
    plt.grid()
    plt.legend()
    plt.xlabel("2011-01-01 Hour")
    plt.ylabel("Aggregated weighted PV value for suitable area")
```
![png](https://github.com/east-winds/geodata/blob/master/images/mask_on_cutout_workflow/output_47_0.png)

```python

```
