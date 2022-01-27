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

In this tutorial, we use the `china` mask, created in this documentation: [mask_creation_workflow](https://github.com/east-winds/geodata/blob/mask/doc/mask/mask_creation_workflow.md)

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
![png](https://github.com/east-winds/geodata/blob/mask/images/mask_on_cutout_workflow/output_21_1.png)

### 4.2 Adding area variable

To calculate and add the variation of grid cell areas by latitude to the cutout, use the `add_grid_area` method. Keeping track of the area for each grid cell is necessary for analyses such as calculating the weighted sum of the grid cells based on their area.

```python
cutout.add_grid_area()
```
### 4.3 Creating PV data through cutout conversion

The code block below will use the `geodata.convert.pv` method to generate `ds_cutout`, an xarray Dataset that contains the pv variable for the cutout. We will also reset the variables names such as `x` and `y` to `lon` and `lat`.

```python
ds_solar = geodata.convert.pv(cutout, panel = "KANENA", orientation = "latitude_optimal")
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
<div><svg style="position: absolute; width: 0; height: 0; overflow: hidden">
<defs>
<symbol id="icon-database" viewBox="0 0 32 32">
<path d="M16 0c-8.837 0-16 2.239-16 5v4c0 2.761 7.163 5 16 5s16-2.239 16-5v-4c0-2.761-7.163-5-16-5z"></path>
<path d="M16 17c-8.837 0-16-2.239-16-5v6c0 2.761 7.163 5 16 5s16-2.239 16-5v-6c0 2.761-7.163 5-16 5z"></path>
<path d="M16 26c-8.837 0-16-2.239-16-5v6c0 2.761 7.163 5 16 5s16-2.239 16-5v-6c0 2.761-7.163 5-16 5z"></path>
</symbol>
<symbol id="icon-file-text2" viewBox="0 0 32 32">
<path d="M28.681 7.159c-0.694-0.947-1.662-2.053-2.724-3.116s-2.169-2.030-3.116-2.724c-1.612-1.182-2.393-1.319-2.841-1.319h-15.5c-1.378 0-2.5 1.121-2.5 2.5v27c0 1.378 1.122 2.5 2.5 2.5h23c1.378 0 2.5-1.122 2.5-2.5v-19.5c0-0.448-0.137-1.23-1.319-2.841zM24.543 5.457c0.959 0.959 1.712 1.825 2.268 2.543h-4.811v-4.811c0.718 0.556 1.584 1.309 2.543 2.268zM28 29.5c0 0.271-0.229 0.5-0.5 0.5h-23c-0.271 0-0.5-0.229-0.5-0.5v-27c0-0.271 0.229-0.5 0.5-0.5 0 0 15.499-0 15.5 0v7c0 0.552 0.448 1 1 1h7v19.5z"></path>
<path d="M23 26h-14c-0.552 0-1-0.448-1-1s0.448-1 1-1h14c0.552 0 1 0.448 1 1s-0.448 1-1 1z"></path>
<path d="M23 22h-14c-0.552 0-1-0.448-1-1s0.448-1 1-1h14c0.552 0 1 0.448 1 1s-0.448 1-1 1z"></path>
<path d="M23 18h-14c-0.552 0-1-0.448-1-1s0.448-1 1-1h14c0.552 0 1 0.448 1 1s-0.448 1-1 1z"></path>
</symbol>
</defs>
</svg>
<style>/* CSS stylesheet for displaying xarray objects in jupyterlab.
 *
 */

:root {
--xr-font-color0: var(--jp-content-font-color0, rgba(0, 0, 0, 1));
--xr-font-color2: var(--jp-content-font-color2, rgba(0, 0, 0, 0.54));
--xr-font-color3: var(--jp-content-font-color3, rgba(0, 0, 0, 0.38));
--xr-border-color: var(--jp-border-color2, #e0e0e0);
--xr-disabled-color: var(--jp-layout-color3, #bdbdbd);
--xr-background-color: var(--jp-layout-color0, white);
--xr-background-color-row-even: var(--jp-layout-color1, white);
--xr-background-color-row-odd: var(--jp-layout-color2, #eeeeee);
}

html[theme=dark],
body.vscode-dark {
--xr-font-color0: rgba(255, 255, 255, 1);
--xr-font-color2: rgba(255, 255, 255, 0.54);
--xr-font-color3: rgba(255, 255, 255, 0.38);
--xr-border-color: #1F1F1F;
--xr-disabled-color: #515151;
--xr-background-color: #111111;
--xr-background-color-row-even: #111111;
--xr-background-color-row-odd: #313131;
}

.xr-wrap {
display: block;
min-width: 300px;
max-width: 700px;
}

.xr-text-repr-fallback {
/* fallback to plain text repr when CSS is not injected (untrusted notebook) */
display: none;
}

.xr-header {
padding-top: 6px;
padding-bottom: 6px;
margin-bottom: 4px;
border-bottom: solid 1px var(--xr-border-color);
}

.xr-header > div,
.xr-header > ul {
display: inline;
margin-top: 0;
margin-bottom: 0;
}

.xr-obj-type,
.xr-array-name {
margin-left: 2px;
margin-right: 10px;
}

.xr-obj-type {
color: var(--xr-font-color2);
}

.xr-sections {
padding-left: 0 !important;
display: grid;
grid-template-columns: 150px auto auto 1fr 20px 20px;
}

.xr-section-item {
display: contents;
}

.xr-section-item input {
display: none;
}

.xr-section-item input + label {
color: var(--xr-disabled-color);
}

.xr-section-item input:enabled + label {
cursor: pointer;
color: var(--xr-font-color2);
}

.xr-section-item input:enabled + label:hover {
color: var(--xr-font-color0);
}

.xr-section-summary {
grid-column: 1;
color: var(--xr-font-color2);
font-weight: 500;
}

.xr-section-summary > span {
display: inline-block;
padding-left: 0.5em;
}

.xr-section-summary-in:disabled + label {
color: var(--xr-font-color2);
}

.xr-section-summary-in + label:before {
display: inline-block;
content: '►';
font-size: 11px;
width: 15px;
text-align: center;
}

.xr-section-summary-in:disabled + label:before {
color: var(--xr-disabled-color);
}

.xr-section-summary-in:checked + label:before {
content: '▼';
}

.xr-section-summary-in:checked + label > span {
display: none;
}

.xr-section-summary,
.xr-section-inline-details {
padding-top: 4px;
padding-bottom: 4px;
}

.xr-section-inline-details {
grid-column: 2 / -1;
}

.xr-section-details {
display: none;
grid-column: 1 / -1;
margin-bottom: 5px;
}

.xr-section-summary-in:checked ~ .xr-section-details {
display: contents;
}

.xr-array-wrap {
grid-column: 1 / -1;
display: grid;
grid-template-columns: 20px auto;
}

.xr-array-wrap > label {
grid-column: 1;
vertical-align: top;
}

.xr-preview {
color: var(--xr-font-color3);
}

.xr-array-preview,
.xr-array-data {
padding: 0 5px !important;
grid-column: 2;
}

.xr-array-data,
.xr-array-in:checked ~ .xr-array-preview {
display: none;
}

.xr-array-in:checked ~ .xr-array-data,
.xr-array-preview {
display: inline-block;
}

.xr-dim-list {
display: inline-block !important;
list-style: none;
padding: 0 !important;
margin: 0;
}

.xr-dim-list li {
display: inline-block;
padding: 0;
margin: 0;
}

.xr-dim-list:before {
content: '(';
}

.xr-dim-list:after {
content: ')';
}

.xr-dim-list li:not(:last-child):after {
content: ',';
padding-right: 5px;
}

.xr-has-index {
font-weight: bold;
}

.xr-var-list,
.xr-var-item {
display: contents;
}

.xr-var-item > div,
.xr-var-item label,
.xr-var-item > .xr-var-name span {
background-color: var(--xr-background-color-row-even);
margin-bottom: 0;
}

.xr-var-item > .xr-var-name:hover span {
padding-right: 5px;
}

.xr-var-list > li:nth-child(odd) > div,
.xr-var-list > li:nth-child(odd) > label,
.xr-var-list > li:nth-child(odd) > .xr-var-name span {
background-color: var(--xr-background-color-row-odd);
}

.xr-var-name {
grid-column: 1;
}

.xr-var-dims {
grid-column: 2;
}

.xr-var-dtype {
grid-column: 3;
text-align: right;
color: var(--xr-font-color2);
}

.xr-var-preview {
grid-column: 4;
}

.xr-var-name,
.xr-var-dims,
.xr-var-dtype,
.xr-preview,
.xr-attrs dt {
white-space: nowrap;
overflow: hidden;
text-overflow: ellipsis;
padding-right: 10px;
}

.xr-var-name:hover,
.xr-var-dims:hover,
.xr-var-dtype:hover,
.xr-attrs dt:hover {
overflow: visible;
width: auto;
z-index: 1;
}

.xr-var-attrs,
.xr-var-data {
display: none;
background-color: var(--xr-background-color) !important;
padding-bottom: 5px !important;
}

.xr-var-attrs-in:checked ~ .xr-var-attrs,
.xr-var-data-in:checked ~ .xr-var-data {
display: block;
}

.xr-var-data > table {
float: right;
}

.xr-var-name span,
.xr-var-data,
.xr-attrs {
padding-left: 25px !important;
}

.xr-attrs,
.xr-var-attrs,
.xr-var-data {
grid-column: 1 / -1;
}

dl.xr-attrs {
padding: 0;
margin: 0;
display: grid;
grid-template-columns: 125px auto;
}

.xr-attrs dt, dd {
padding: 0;
margin: 0;
float: left;
padding-right: 10px;
width: auto;
}

.xr-attrs dt {
font-weight: normal;
grid-column: 1;
}

.xr-attrs dt:hover span {
display: inline-block;
background: var(--xr-background-color);
padding-right: 10px;
}

.xr-attrs dd {
grid-column: 2;
white-space: pre-wrap;
word-break: break-all;
}

.xr-icon-database,
.xr-icon-file-text2 {
display: inline-block;
vertical-align: middle;
width: 1em;
height: 1.5em !important;
stroke-width: 0;
stroke: currentColor;
fill: currentColor;
}
</style><pre class='xr-text-repr-fallback'>&lt;xarray.Dataset&gt;
Dimensions:  (lat: 73, lon: 101, time: 1)
Coordinates:

* lat      (lat) float64 18.0 18.5 19.0 19.5 20.0 ... 52.0 52.5 53.0 53.5 54.0
* time     (time) datetime64[ns] 2011-01-01T12:00:00
* lon      (lon) float64 73.12 73.75 74.38 75.0 ... 133.8 134.4 135.0 135.6
  Data variables:
  solar    (lat, time, lon) float64 0.1928 0.1917 0.1641 ... 0.06296 0.04254
  mask     (lat, lon) float32 0.0 0.0 0.0 0.0 0.0 0.0 ... 0.0 0.0 0.0 0.0 0.0
  area     (lat, lon) float64 3.663e+03 3.663e+03 ... 2.281e+03 2.281e+03</pre><div class='xr-wrap' hidden><div class='xr-header'><div class='xr-obj-type'>xarray.Dataset</div></div><ul class='xr-sections'><li class='xr-section-item'><input id='section-8ab24a01-1a74-46f9-a905-7a491de72ce0' class='xr-section-summary-in' type='checkbox' disabled ><label for='section-8ab24a01-1a74-46f9-a905-7a491de72ce0' class='xr-section-summary'  title='Expand/collapse section'>Dimensions:</label><div class='xr-section-inline-details'><ul class='xr-dim-list'><li><span class='xr-has-index'>lat</span>: 73</li><li><span class='xr-has-index'>lon</span>: 101</li><li><span class='xr-has-index'>time</span>: 1</li></ul></div><div class='xr-section-details'></div></li><li class='xr-section-item'><input id='section-9547bf7f-9af7-4f14-8d8a-47b1178823ff' class='xr-section-summary-in' type='checkbox'  checked><label for='section-9547bf7f-9af7-4f14-8d8a-47b1178823ff' class='xr-section-summary' >Coordinates: <span>(3)</span></label><div class='xr-section-inline-details'></div><div class='xr-section-details'><ul class='xr-var-list'><li class='xr-var-item'><div class='xr-var-name'><span class='xr-has-index'>lat</span></div><div class='xr-var-dims'>(lat)</div><div class='xr-var-dtype'>float64</div><div class='xr-var-preview xr-preview'>18.0 18.5 19.0 ... 53.0 53.5 54.0</div><input id='attrs-d463eb31-c61c-406b-a4b8-e323cafb5e49' class='xr-var-attrs-in' type='checkbox' ><label for='attrs-d463eb31-c61c-406b-a4b8-e323cafb5e49' title='Show/Hide attributes'><svg class='icon xr-icon-file-text2'><use xlink:href='#icon-file-text2'></use></svg></label><input id='data-c07f8418-7e04-4bbd-ae5b-7e1dc417127f' class='xr-var-data-in' type='checkbox'><label for='data-c07f8418-7e04-4bbd-ae5b-7e1dc417127f' title='Show/Hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-var-attrs'><dl class='xr-attrs'><dt><span>long_name :</span></dt><dd>latitude</dd><dt><span>units :</span></dt><dd>degrees_north</dd><dt><span>vmax :</span></dt><dd>1000000000000000.0</dd><dt><span>vmin :</span></dt><dd>-1000000000000000.0</dd><dt><span>valid_range :</span></dt><dd>[-1.e+15  1.e+15]</dd></dl></div><div class='xr-var-data'><pre>array([18. , 18.5, 19. , 19.5, 20. , 20.5, 21. , 21.5, 22. , 22.5, 23. , 23.5,
  24. , 24.5, 25. , 25.5, 26. , 26.5, 27. , 27.5, 28. , 28.5, 29. , 29.5,
  30. , 30.5, 31. , 31.5, 32. , 32.5, 33. , 33.5, 34. , 34.5, 35. , 35.5,
  36. , 36.5, 37. , 37.5, 38. , 38.5, 39. , 39.5, 40. , 40.5, 41. , 41.5,
  42. , 42.5, 43. , 43.5, 44. , 44.5, 45. , 45.5, 46. , 46.5, 47. , 47.5,
  48. , 48.5, 49. , 49.5, 50. , 50.5, 51. , 51.5, 52. , 52.5, 53. , 53.5,
  54. ])</pre></div></li><li class='xr-var-item'><div class='xr-var-name'><span class='xr-has-index'>time</span></div><div class='xr-var-dims'>(time)</div><div class='xr-var-dtype'>datetime64[ns]</div><div class='xr-var-preview xr-preview'>2011-01-01T12:00:00</div><input id='attrs-95496c56-82d4-46cb-85f2-fcdf31790bf1' class='xr-var-attrs-in' type='checkbox' disabled><label for='attrs-95496c56-82d4-46cb-85f2-fcdf31790bf1' title='Show/Hide attributes'><svg class='icon xr-icon-file-text2'><use xlink:href='#icon-file-text2'></use></svg></label><input id='data-ab4e17b6-3ac8-4ca5-9117-d1f2553f37f9' class='xr-var-data-in' type='checkbox'><label for='data-ab4e17b6-3ac8-4ca5-9117-d1f2553f37f9' title='Show/Hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-var-attrs'><dl class='xr-attrs'></dl></div><div class='xr-var-data'><pre>array([&#x27;2011-01-01T12:00:00.000000000&#x27;], dtype=&#x27;datetime64[ns]&#x27;)</pre></div></li><li class='xr-var-item'><div class='xr-var-name'><span class='xr-has-index'>lon</span></div><div class='xr-var-dims'>(lon)</div><div class='xr-var-dtype'>float64</div><div class='xr-var-preview xr-preview'>73.12 73.75 74.38 ... 135.0 135.6</div><input id='attrs-47190f7b-ce3e-4426-ae66-6e970f6aa08d' class='xr-var-attrs-in' type='checkbox' ><label for='attrs-47190f7b-ce3e-4426-ae66-6e970f6aa08d' title='Show/Hide attributes'><svg class='icon xr-icon-file-text2'><use xlink:href='#icon-file-text2'></use></svg></label><input id='data-b1239ae7-1077-40a9-b4db-3334b00cbce1' class='xr-var-data-in' type='checkbox'><label for='data-b1239ae7-1077-40a9-b4db-3334b00cbce1' title='Show/Hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-var-attrs'><dl class='xr-attrs'><dt><span>long_name :</span></dt><dd>longitude</dd><dt><span>units :</span></dt><dd>degrees_east</dd><dt><span>vmax :</span></dt><dd>1000000000000000.0</dd><dt><span>vmin :</span></dt><dd>-1000000000000000.0</dd><dt><span>valid_range :</span></dt><dd>[-1.e+15  1.e+15]</dd></dl></div><div class='xr-var-data'><pre>array([ 73.125,  73.75 ,  74.375,  75.   ,  75.625,  76.25 ,  76.875,  77.5  ,
  78.125,  78.75 ,  79.375,  80.   ,  80.625,  81.25 ,  81.875,  82.5  ,
  83.125,  83.75 ,  84.375,  85.   ,  85.625,  86.25 ,  86.875,  87.5  ,
  88.125,  88.75 ,  89.375,  90.   ,  90.625,  91.25 ,  91.875,  92.5  ,
  93.125,  93.75 ,  94.375,  95.   ,  95.625,  96.25 ,  96.875,  97.5  ,
  98.125,  98.75 ,  99.375, 100.   , 100.625, 101.25 , 101.875, 102.5  ,
  103.125, 103.75 , 104.375, 105.   , 105.625, 106.25 , 106.875, 107.5  ,
  108.125, 108.75 , 109.375, 110.   , 110.625, 111.25 , 111.875, 112.5  ,
  113.125, 113.75 , 114.375, 115.   , 115.625, 116.25 , 116.875, 117.5  ,
  118.125, 118.75 , 119.375, 120.   , 120.625, 121.25 , 121.875, 122.5  ,
  123.125, 123.75 , 124.375, 125.   , 125.625, 126.25 , 126.875, 127.5  ,
  128.125, 128.75 , 129.375, 130.   , 130.625, 131.25 , 131.875, 132.5  ,
  133.125, 133.75 , 134.375, 135.   , 135.625])</pre></div></li></ul></div></li><li class='xr-section-item'><input id='section-0382ca35-f497-4554-9576-507a54960721' class='xr-section-summary-in' type='checkbox'  checked><label for='section-0382ca35-f497-4554-9576-507a54960721' class='xr-section-summary' >Data variables: <span>(3)</span></label><div class='xr-section-inline-details'></div><div class='xr-section-details'><ul class='xr-var-list'><li class='xr-var-item'><div class='xr-var-name'><span>solar</span></div><div class='xr-var-dims'>(lat, time, lon)</div><div class='xr-var-dtype'>float64</div><div class='xr-var-preview xr-preview'>0.1928 0.1917 ... 0.06296 0.04254</div><input id='attrs-aed474e2-277d-40fa-ad89-69d2796e7f73' class='xr-var-attrs-in' type='checkbox' disabled><label for='attrs-aed474e2-277d-40fa-ad89-69d2796e7f73' title='Show/Hide attributes'><svg class='icon xr-icon-file-text2'><use xlink:href='#icon-file-text2'></use></svg></label><input id='data-7e58d5d5-1675-460e-96e3-ccfb5c4c4782' class='xr-var-data-in' type='checkbox'><label for='data-7e58d5d5-1675-460e-96e3-ccfb5c4c4782' title='Show/Hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-var-attrs'><dl class='xr-attrs'></dl></div><div class='xr-var-data'><pre>array([[[0.19283223, 0.19170759, 0.16413559, ..., 0.10883634,
  0.1237414 , 0.14484597]],

  [[0.18934099, 0.17642454, 0.15045152, ..., 0.11743452,
  0.12400376, 0.13160104]],

  [[0.18668941, 0.16858997, 0.15180354, ..., 0.11671159,
  0.12352632, 0.12845966]],

  ...,

  [[0.08910529, 0.08343294, 0.0793837 , ..., 0.08158544,
  0.07718862, 0.05804729]],

  [[0.08211417, 0.07891589, 0.07478246, ..., 0.06911246,
  0.08452805, 0.07724436]],

  [[0.07631315, 0.07576691, 0.07252836, ..., 0.08120154,
  0.06295666, 0.04253606]]])</pre></div></li><li class='xr-var-item'><div class='xr-var-name'><span>mask</span></div><div class='xr-var-dims'>(lat, lon)</div><div class='xr-var-dtype'>float32</div><div class='xr-var-preview xr-preview'>0.0 0.0 0.0 0.0 ... 0.0 0.0 0.0 0.0</div><input id='attrs-349f5593-4b9f-4c98-896f-826bb6d64e4e' class='xr-var-attrs-in' type='checkbox' disabled><label for='attrs-349f5593-4b9f-4c98-896f-826bb6d64e4e' title='Show/Hide attributes'><svg class='icon xr-icon-file-text2'><use xlink:href='#icon-file-text2'></use></svg></label><input id='data-c6b1b162-bed0-4c9a-a49a-e1cada0d4624' class='xr-var-data-in' type='checkbox'><label for='data-c6b1b162-bed0-4c9a-a49a-e1cada0d4624' title='Show/Hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-var-attrs'><dl class='xr-attrs'></dl></div><div class='xr-var-data'><pre>array([[0., 0., 0., ..., 0., 0., 0.],
  [0., 0., 0., ..., 0., 0., 0.],
  [0., 0., 0., ..., 0., 0., 0.],
  ...,
  [0., 0., 0., ..., 0., 0., 0.],
  [0., 0., 0., ..., 0., 0., 0.],
  [0., 0., 0., ..., 0., 0., 0.]], dtype=float32)</pre></div></li><li class='xr-var-item'><div class='xr-var-name'><span>area</span></div><div class='xr-var-dims'>(lat, lon)</div><div class='xr-var-dtype'>float64</div><div class='xr-var-preview xr-preview'>3.663e+03 3.663e+03 ... 2.281e+03</div><input id='attrs-ff072cfd-7ef5-4621-bd42-512fc973a656' class='xr-var-attrs-in' type='checkbox' disabled><label for='attrs-ff072cfd-7ef5-4621-bd42-512fc973a656' title='Show/Hide attributes'><svg class='icon xr-icon-file-text2'><use xlink:href='#icon-file-text2'></use></svg></label><input id='data-5137f757-331e-4361-ac41-39a5a9a45b9c' class='xr-var-data-in' type='checkbox'><label for='data-5137f757-331e-4361-ac41-39a5a9a45b9c' title='Show/Hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-var-attrs'><dl class='xr-attrs'></dl></div><div class='xr-var-data'><pre>array([[3662.98, 3662.98, 3662.98, ..., 3662.98, 3662.98, 3662.98],
  [3652.71, 3652.71, 3652.71, ..., 3652.71, 3652.71, 3652.71],
  [3642.17, 3642.17, 3642.17, ..., 3642.17, 3642.17, 3642.17],
  ...,
  [2334.79, 2334.79, 2334.79, ..., 2334.79, 2334.79, 2334.79],
  [2307.92, 2307.92, 2307.92, ..., 2307.92, 2307.92, 2307.92],
  [2280.87, 2280.87, 2280.87, ..., 2280.87, 2280.87, 2280.87]])</pre></div></li></ul></div></li><li class='xr-section-item'><input id='section-d9ae16b3-ba53-403a-9a76-1f03e90915ad' class='xr-section-summary-in' type='checkbox' disabled ><label for='section-d9ae16b3-ba53-403a-9a76-1f03e90915ad' class='xr-section-summary'  title='Expand/collapse section'>Attributes: <span>(0)</span></label><div class='xr-section-inline-details'></div><div class='xr-section-details'><dl class='xr-attrs'></dl></div></li></ul></div></div>

Visualize the averaged PV value for each grid cell in the Cutout. Note that the data is the aggregated value for the date.

```python
combine_mean['Jiangsu']['solar'].plot()
```
<matplotlib.collections.QuadMesh at 0x1ee15dd8a90>
![png](https://github.com/east-winds/geodata/blob/mask/images/mask_on_cutout_workflow/output_37_1.png)

Visualize the masking value for each grid cell in the Cutout.

```python
combine_mean['Jiangsu']['mask'].plot()
```
<matplotlib.collections.QuadMesh at 0x1ee15e72520>
![png](https://github.com/east-winds/geodata/blob/mask/images/mask_on_cutout_workflow/output_39_1.png)

**Area- and mask-weighted hourly PV values**

We use the raw hourly output generated by cutout to create time-series PV plots weighted by the mask and area. Note that we transposed ds_cutout so that time is set as the first dimension, which ease the following calculation since we want to aggregate the array spatially from each grid cell.

```python
combine = cutout.mask(dataset = ds_cutout.transpose("time", "lat", "lon"))
combine.keys()
```
```dict_keys(['merged_mask', 'Jiangsu', 'Shanghai', 'Zhejiang'])```


Calculate the aggregated mean solar PV for each provinces, at each time point. We will apply this equation below to calculate the area-weighted average. We save the result into a dictionary `PV_dict`, where its keys are the provinces, and the corresponding values are the PV series.

\begin{equation}
Aggregated\:Solar\:Power\:for\:each\:region = \frac{\sum_{}^{for\:each\:grid\:cell}Grid\:cell\:area*mask\:value*Solar\:power}{\sum_{}^{for\:each\:grid\:cell}Grid\:cell\:area*mask\:value}
\end{equation}

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
<div><svg style="position: absolute; width: 0; height: 0; overflow: hidden">
<defs>
<symbol id="icon-database" viewBox="0 0 32 32">
<path d="M16 0c-8.837 0-16 2.239-16 5v4c0 2.761 7.163 5 16 5s16-2.239 16-5v-4c0-2.761-7.163-5-16-5z"></path>
<path d="M16 17c-8.837 0-16-2.239-16-5v6c0 2.761 7.163 5 16 5s16-2.239 16-5v-6c0 2.761-7.163 5-16 5z"></path>
<path d="M16 26c-8.837 0-16-2.239-16-5v6c0 2.761 7.163 5 16 5s16-2.239 16-5v-6c0 2.761-7.163 5-16 5z"></path>
</symbol>
<symbol id="icon-file-text2" viewBox="0 0 32 32">
<path d="M28.681 7.159c-0.694-0.947-1.662-2.053-2.724-3.116s-2.169-2.030-3.116-2.724c-1.612-1.182-2.393-1.319-2.841-1.319h-15.5c-1.378 0-2.5 1.121-2.5 2.5v27c0 1.378 1.122 2.5 2.5 2.5h23c1.378 0 2.5-1.122 2.5-2.5v-19.5c0-0.448-0.137-1.23-1.319-2.841zM24.543 5.457c0.959 0.959 1.712 1.825 2.268 2.543h-4.811v-4.811c0.718 0.556 1.584 1.309 2.543 2.268zM28 29.5c0 0.271-0.229 0.5-0.5 0.5h-23c-0.271 0-0.5-0.229-0.5-0.5v-27c0-0.271 0.229-0.5 0.5-0.5 0 0 15.499-0 15.5 0v7c0 0.552 0.448 1 1 1h7v19.5z"></path>
<path d="M23 26h-14c-0.552 0-1-0.448-1-1s0.448-1 1-1h14c0.552 0 1 0.448 1 1s-0.448 1-1 1z"></path>
<path d="M23 22h-14c-0.552 0-1-0.448-1-1s0.448-1 1-1h14c0.552 0 1 0.448 1 1s-0.448 1-1 1z"></path>
<path d="M23 18h-14c-0.552 0-1-0.448-1-1s0.448-1 1-1h14c0.552 0 1 0.448 1 1s-0.448 1-1 1z"></path>
</symbol>
</defs>
</svg>
<style>/* CSS stylesheet for displaying xarray objects in jupyterlab.
 *
 */

:root {
--xr-font-color0: var(--jp-content-font-color0, rgba(0, 0, 0, 1));
--xr-font-color2: var(--jp-content-font-color2, rgba(0, 0, 0, 0.54));
--xr-font-color3: var(--jp-content-font-color3, rgba(0, 0, 0, 0.38));
--xr-border-color: var(--jp-border-color2, #e0e0e0);
--xr-disabled-color: var(--jp-layout-color3, #bdbdbd);
--xr-background-color: var(--jp-layout-color0, white);
--xr-background-color-row-even: var(--jp-layout-color1, white);
--xr-background-color-row-odd: var(--jp-layout-color2, #eeeeee);
}

html[theme=dark],
body.vscode-dark {
--xr-font-color0: rgba(255, 255, 255, 1);
--xr-font-color2: rgba(255, 255, 255, 0.54);
--xr-font-color3: rgba(255, 255, 255, 0.38);
--xr-border-color: #1F1F1F;
--xr-disabled-color: #515151;
--xr-background-color: #111111;
--xr-background-color-row-even: #111111;
--xr-background-color-row-odd: #313131;
}

.xr-wrap {
display: block;
min-width: 300px;
max-width: 700px;
}

.xr-text-repr-fallback {
/* fallback to plain text repr when CSS is not injected (untrusted notebook) */
display: none;
}

.xr-header {
padding-top: 6px;
padding-bottom: 6px;
margin-bottom: 4px;
border-bottom: solid 1px var(--xr-border-color);
}

.xr-header > div,
.xr-header > ul {
display: inline;
margin-top: 0;
margin-bottom: 0;
}

.xr-obj-type,
.xr-array-name {
margin-left: 2px;
margin-right: 10px;
}

.xr-obj-type {
color: var(--xr-font-color2);
}

.xr-sections {
padding-left: 0 !important;
display: grid;
grid-template-columns: 150px auto auto 1fr 20px 20px;
}

.xr-section-item {
display: contents;
}

.xr-section-item input {
display: none;
}

.xr-section-item input + label {
color: var(--xr-disabled-color);
}

.xr-section-item input:enabled + label {
cursor: pointer;
color: var(--xr-font-color2);
}

.xr-section-item input:enabled + label:hover {
color: var(--xr-font-color0);
}

.xr-section-summary {
grid-column: 1;
color: var(--xr-font-color2);
font-weight: 500;
}

.xr-section-summary > span {
display: inline-block;
padding-left: 0.5em;
}

.xr-section-summary-in:disabled + label {
color: var(--xr-font-color2);
}

.xr-section-summary-in + label:before {
display: inline-block;
content: '►';
font-size: 11px;
width: 15px;
text-align: center;
}

.xr-section-summary-in:disabled + label:before {
color: var(--xr-disabled-color);
}

.xr-section-summary-in:checked + label:before {
content: '▼';
}

.xr-section-summary-in:checked + label > span {
display: none;
}

.xr-section-summary,
.xr-section-inline-details {
padding-top: 4px;
padding-bottom: 4px;
}

.xr-section-inline-details {
grid-column: 2 / -1;
}

.xr-section-details {
display: none;
grid-column: 1 / -1;
margin-bottom: 5px;
}

.xr-section-summary-in:checked ~ .xr-section-details {
display: contents;
}

.xr-array-wrap {
grid-column: 1 / -1;
display: grid;
grid-template-columns: 20px auto;
}

.xr-array-wrap > label {
grid-column: 1;
vertical-align: top;
}

.xr-preview {
color: var(--xr-font-color3);
}

.xr-array-preview,
.xr-array-data {
padding: 0 5px !important;
grid-column: 2;
}

.xr-array-data,
.xr-array-in:checked ~ .xr-array-preview {
display: none;
}

.xr-array-in:checked ~ .xr-array-data,
.xr-array-preview {
display: inline-block;
}

.xr-dim-list {
display: inline-block !important;
list-style: none;
padding: 0 !important;
margin: 0;
}

.xr-dim-list li {
display: inline-block;
padding: 0;
margin: 0;
}

.xr-dim-list:before {
content: '(';
}

.xr-dim-list:after {
content: ')';
}

.xr-dim-list li:not(:last-child):after {
content: ',';
padding-right: 5px;
}

.xr-has-index {
font-weight: bold;
}

.xr-var-list,
.xr-var-item {
display: contents;
}

.xr-var-item > div,
.xr-var-item label,
.xr-var-item > .xr-var-name span {
background-color: var(--xr-background-color-row-even);
margin-bottom: 0;
}

.xr-var-item > .xr-var-name:hover span {
padding-right: 5px;
}

.xr-var-list > li:nth-child(odd) > div,
.xr-var-list > li:nth-child(odd) > label,
.xr-var-list > li:nth-child(odd) > .xr-var-name span {
background-color: var(--xr-background-color-row-odd);
}

.xr-var-name {
grid-column: 1;
}

.xr-var-dims {
grid-column: 2;
}

.xr-var-dtype {
grid-column: 3;
text-align: right;
color: var(--xr-font-color2);
}

.xr-var-preview {
grid-column: 4;
}

.xr-var-name,
.xr-var-dims,
.xr-var-dtype,
.xr-preview,
.xr-attrs dt {
white-space: nowrap;
overflow: hidden;
text-overflow: ellipsis;
padding-right: 10px;
}

.xr-var-name:hover,
.xr-var-dims:hover,
.xr-var-dtype:hover,
.xr-attrs dt:hover {
overflow: visible;
width: auto;
z-index: 1;
}

.xr-var-attrs,
.xr-var-data {
display: none;
background-color: var(--xr-background-color) !important;
padding-bottom: 5px !important;
}

.xr-var-attrs-in:checked ~ .xr-var-attrs,
.xr-var-data-in:checked ~ .xr-var-data {
display: block;
}

.xr-var-data > table {
float: right;
}

.xr-var-name span,
.xr-var-data,
.xr-attrs {
padding-left: 25px !important;
}

.xr-attrs,
.xr-var-attrs,
.xr-var-data {
grid-column: 1 / -1;
}

dl.xr-attrs {
padding: 0;
margin: 0;
display: grid;
grid-template-columns: 125px auto;
}

.xr-attrs dt, dd {
padding: 0;
margin: 0;
float: left;
padding-right: 10px;
width: auto;
}

.xr-attrs dt {
font-weight: normal;
grid-column: 1;
}

.xr-attrs dt:hover span {
display: inline-block;
background: var(--xr-background-color);
padding-right: 10px;
}

.xr-attrs dd {
grid-column: 2;
white-space: pre-wrap;
word-break: break-all;
}

.xr-icon-database,
.xr-icon-file-text2 {
display: inline-block;
vertical-align: middle;
width: 1em;
height: 1.5em !important;
stroke-width: 0;
stroke: currentColor;
fill: currentColor;
}
</style><pre class='xr-text-repr-fallback'>&lt;xarray.DataArray (time: 24)&gt;
array([0.17420805, 0.36234572, 0.51907652, 0.62608188, 0.67846743,
0.66500807, 0.56758925, 0.38363103, 0.1380775 , 0.        ,
0.        , 0.        , 0.        , 0.        , 0.        ,
0.        , 0.        , 0.        , 0.        , 0.        ,
0.        , 0.        , 0.        , 0.03283585])
Coordinates:

* time     (time) datetime64[ns] 2011-01-01T00:30:00 ... 2011-01-01T23:30:00</pre><div class='xr-wrap' hidden><div class='xr-header'><div class='xr-obj-type'>xarray.DataArray</div><div class='xr-array-name'></div><ul class='xr-dim-list'><li><span class='xr-has-index'>time</span>: 24</li></ul></div><ul class='xr-sections'><li class='xr-section-item'><div class='xr-array-wrap'><input id='section-b5857dc0-baf1-4491-b908-651b954e49ba' class='xr-array-in' type='checkbox' checked><label for='section-b5857dc0-baf1-4491-b908-651b954e49ba' title='Show/hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-array-preview xr-preview'><span>0.1742 0.3623 0.5191 0.6261 0.6785 0.665 ... 0.0 0.0 0.0 0.0 0.03284</span></div><div class='xr-array-data'><pre>array([0.17420805, 0.36234572, 0.51907652, 0.62608188, 0.67846743,
  0.66500807, 0.56758925, 0.38363103, 0.1380775 , 0.        ,
  0.        , 0.        , 0.        , 0.        , 0.        ,
  0.        , 0.        , 0.        , 0.        , 0.        ,
  0.        , 0.        , 0.        , 0.03283585])</pre></div></div></li><li class='xr-section-item'><input id='section-9ab591cc-1319-40ad-8c5a-97b0df367b3b' class='xr-section-summary-in' type='checkbox'  checked><label for='section-9ab591cc-1319-40ad-8c5a-97b0df367b3b' class='xr-section-summary' >Coordinates: <span>(1)</span></label><div class='xr-section-inline-details'></div><div class='xr-section-details'><ul class='xr-var-list'><li class='xr-var-item'><div class='xr-var-name'><span class='xr-has-index'>time</span></div><div class='xr-var-dims'>(time)</div><div class='xr-var-dtype'>datetime64[ns]</div><div class='xr-var-preview xr-preview'>2011-01-01T00:30:00 ... 2011-01-...</div><input id='attrs-52a8ca48-9e02-4e11-89f5-4303d8957494' class='xr-var-attrs-in' type='checkbox' disabled><label for='attrs-52a8ca48-9e02-4e11-89f5-4303d8957494' title='Show/Hide attributes'><svg class='icon xr-icon-file-text2'><use xlink:href='#icon-file-text2'></use></svg></label><input id='data-11be99ba-6e70-49b6-9dc9-9b79a5f7f4ff' class='xr-var-data-in' type='checkbox'><label for='data-11be99ba-6e70-49b6-9dc9-9b79a5f7f4ff' title='Show/Hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-var-attrs'><dl class='xr-attrs'></dl></div><div class='xr-var-data'><pre>array([&#x27;2011-01-01T00:30:00.000000000&#x27;, &#x27;2011-01-01T01:30:00.000000000&#x27;,
  &#x27;2011-01-01T02:30:00.000000000&#x27;, &#x27;2011-01-01T03:30:00.000000000&#x27;,
  &#x27;2011-01-01T04:30:00.000000000&#x27;, &#x27;2011-01-01T05:30:00.000000000&#x27;,
  &#x27;2011-01-01T06:30:00.000000000&#x27;, &#x27;2011-01-01T07:30:00.000000000&#x27;,
  &#x27;2011-01-01T08:30:00.000000000&#x27;, &#x27;2011-01-01T09:30:00.000000000&#x27;,
  &#x27;2011-01-01T10:30:00.000000000&#x27;, &#x27;2011-01-01T11:30:00.000000000&#x27;,
  &#x27;2011-01-01T12:30:00.000000000&#x27;, &#x27;2011-01-01T13:30:00.000000000&#x27;,
  &#x27;2011-01-01T14:30:00.000000000&#x27;, &#x27;2011-01-01T15:30:00.000000000&#x27;,
  &#x27;2011-01-01T16:30:00.000000000&#x27;, &#x27;2011-01-01T17:30:00.000000000&#x27;,
  &#x27;2011-01-01T18:30:00.000000000&#x27;, &#x27;2011-01-01T19:30:00.000000000&#x27;,
  &#x27;2011-01-01T20:30:00.000000000&#x27;, &#x27;2011-01-01T21:30:00.000000000&#x27;,
  &#x27;2011-01-01T22:30:00.000000000&#x27;, &#x27;2011-01-01T23:30:00.000000000&#x27;],
  dtype=&#x27;datetime64[ns]&#x27;)</pre></div></li></ul></div></li><li class='xr-section-item'><input id='section-8292e543-f805-47db-b451-d94bde7b5220' class='xr-section-summary-in' type='checkbox' disabled ><label for='section-8292e543-f805-47db-b451-d94bde7b5220' class='xr-section-summary'  title='Expand/collapse section'>Attributes: <span>(0)</span></label><div class='xr-section-inline-details'></div><div class='xr-section-details'><dl class='xr-attrs'></dl></div></li></ul></div></div>

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
![png](https://github.com/east-winds/geodata/blob/mask/images/mask_on_cutout_workflow/output_47_0.png)

```python

```
