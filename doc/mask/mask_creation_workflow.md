# Mask Workflow Example

## 1. Introduction

Geodata is able to process geospatial data to extract cutouts over specified geographies. Built off the [rasterio library](https://rasterio.readthedocs.io/en/latest/quickstart.html), the **mask** module imports rasters and shapefiles, merges and flattens multiple layers together, and extracts subsetted cutout data from merged masks and shapefiles.

Functionalities explored in this notebook:

- [Creating a mask object, adding and manipulating layers](#3-Creating-mask-object-adding-and-manipulating-layers)
- [Opening a shapefile and adding shape features as layers](#4-Openning-shapefile-and-adding-shape-features-as-a-layer)
- [Merging and flattening layers](#5-Merging-and-flattening-layers)
- [Eliminate small contiguous areas](6-Eliminate-small-contiguous-areas)
- [Extracting shapes from mask](#7-Extracting-shapes-from-mask)
- [Saving and loading masks](#8-Saving-and-Loading-masks)

## 2. Setup

To start, import the geodata package and required libraries. We can also import the `geodata.mask.show()` method for simplicity of its use.

```python
import geodata
import numpy as np
from geodata.mask import show 
import matplotlib.pyplot as plt
import geopandas as gpd
```

To launch a logger for detailed debugging, run:

```python
import logging
logging.basicConfig(level=logging.INFO)
```

We use [cartopy](https://scitools.org.uk/cartopy/docs/latest/tutorials/using_the_shapereader.html#cartopy.io.shapereader.Reader) to download some common administrative region shapes, but user-provided shapefiles will also work:

```python
import cartopy.io.shapereader as shpreader
```

We will use the following geotiff and shape files for this demo:

#### china_modis.tif

We downloaded the MODIS land cover data, which uses satellite remote sensing data to estimate the land use type on an annual basis. See: [EarthData_MCD12Q1](https://lpdaac.usgs.gov/products/mcd12q1v006/).

We will use the IGBP classification ('LC_Type1') which has 17 different land use characterizations (the corresponding data thus takes values from 1.0 to 17.0).
All the "Bands" are listed here: [Google_earth_engine_MODIS_006_MCD12Q1](https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MCD12Q1#bands)

#### china_elevation.tif & china_slope.tif

These two rasters are based on the elevation map from: [Google_earth_engine_MODIS_CGIAR_SRTM90_V4](https://developers.google.com/earth-engine/datasets/catalog/CGIAR_SRTM90_V4?hl=en). Slope was computed in degrees using the 4-connected neighbors of each pixel.

#### UNEP_WDPA_China shapefiles

We downloaded the environmental protected area for China from: [ProtectedPlanet_China](https://www.protectedplanet.net/country/CHN). These shapefiles are distributed among 3 subfolders upon successful download and decompression due to the large size. We will create path variables for all three subfolders and we will only take the polygon shapes.

Alternatively, We can also retrieve the environmental protected area from Google Earth Engine: [Google_earth_engine_WCMC_WDPA](https://developers.google.com/earth-engine/datasets/catalog/WCMC_WDPA_current_polygons). The shapefile will contain the protected shapes from entire world (and the size is slightly over 1 GB), and additional data cleaning will be necessary if the user wants just the shapes within China.

```python
modis_path = 'data/china_modis.tif'
elevation_path = 'data/china_elevation.tif'
slope_path = 'data/china_slope.tif'

wdpa_shape_path_0 = 'data/shapefiles/0/WDPA_WDOECM_Nov2021_Public_CHN_shp-polygons.shp'
wdpa_shape_path_1 = 'data/shapefiles/1/WDPA_WDOECM_Nov2021_Public_CHN_shp-polygons.shp'
wdpa_shape_path_2 = 'data/shapefiles/2/WDPA_WDOECM_Nov2021_Public_CHN_shp-polygons.shp'
```

Let us get province shapes from `cartopy` and save the path as `prov_path`. This can also be the path to user-supplied shape files.

```python
prov_path = shpreader.natural_earth(resolution='10m', category='cultural', 
                                    name = 'admin_1_states_provinces')
prov_path
```

'C:\\Users\\fengj\\.local\\share\\cartopy\\shapefiles\\natural_earth\\cultural\\ne_10m_admin_1_states_provinces.shp'
Load the shapes contained in path `prov_path` using the `geopandas` library.

```python
all_shapes = gpd.read_file(prov_path, encoding = 'utf-8')
all_shapes.head(2)
```
<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>featurecla</th>
      <th>scalerank</th>
      <th>adm1_code</th>
      <th>diss_me</th>
      <th>iso_3166_2</th>
      <th>wikipedia</th>
      <th>iso_a2</th>
      <th>adm0_sr</th>
      <th>name</th>
      <th>name_alt</th>
      <th>...</th>
      <th>name_nl</th>
      <th>name_pl</th>
      <th>name_pt</th>
      <th>name_ru</th>
      <th>name_sv</th>
      <th>name_tr</th>
      <th>name_vi</th>
      <th>name_zh</th>
      <th>ne_id</th>
      <th>geometry</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>Admin-1 scale rank</td>
      <td>3</td>
      <td>ARG-1309</td>
      <td>1309</td>
      <td>AR-E</td>
      <td>None</td>
      <td>AR</td>
      <td>1</td>
      <td>Entre Ríos</td>
      <td>Entre-Rios</td>
      <td>...</td>
      <td>Entre Ríos</td>
      <td>Entre Ríos</td>
      <td>Entre Ríos</td>
      <td>Энтре-Риос</td>
      <td>Entre Ríos</td>
      <td>Entre Ríos eyaleti</td>
      <td>Entre Ríos</td>
      <td>恩特雷里奥斯省</td>
      <td>1159309789</td>
      <td>POLYGON ((-58.20011 -32.44713, -58.20012 -32.4...</td>
    </tr>
    <tr>
      <th>1</th>
      <td>Admin-1 scale rank</td>
      <td>6</td>
      <td>URY-8</td>
      <td>8</td>
      <td>UY-PA</td>
      <td>None</td>
      <td>UY</td>
      <td>1</td>
      <td>Paysandú</td>
      <td>None</td>
      <td>...</td>
      <td>Paysandú</td>
      <td>Paysandú</td>
      <td>Paysandú</td>
      <td>Пайсанду</td>
      <td>Paysandú</td>
      <td>Paysandu Departmanı</td>
      <td>Paysandú</td>
      <td>派桑杜省</td>
      <td>1159307733</td>
      <td>POLYGON ((-58.20012 -32.44720, -58.20011 -32.4...</td>
    </tr>
  </tbody>
</table>
<p>2 rows × 84 columns</p>
</div>

Geopandas data filtering with geodataframe is identical to pandas. Let us select all the rows that contains a Chinese shape.

```python
china_shapes = all_shapes[all_shapes['admin'] == 'China']
```
Next, to load the WDPA environmental protected shapefiles as a layer in the china mask, we will use the GeoPandas library. `gpd.read_file()` will return a GeoPandas dataframe including shape attributes and geometry given the file path. Like Pandas, we can read multiple dataframes and append them together. In the code below, we will create one GeoPandas dataframe from three paths that we have for the Chinese environmental protected shapes.

```python
wdpa_shapes = gpd.read_file(wdpa_shape_path_0
             ).append(gpd.read_file(wdpa_shape_path_1)
             ).append(gpd.read_file(wdpa_shape_path_2))
wdpa_shapes.head(2)
```
<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>WDPAID</th>
      <th>WDPA_PID</th>
      <th>PA_DEF</th>
      <th>NAME</th>
      <th>ORIG_NAME</th>
      <th>DESIG</th>
      <th>DESIG_ENG</th>
      <th>DESIG_TYPE</th>
      <th>IUCN_CAT</th>
      <th>INT_CRIT</th>
      <th>...</th>
      <th>MANG_AUTH</th>
      <th>MANG_PLAN</th>
      <th>VERIF</th>
      <th>METADATAID</th>
      <th>SUB_LOC</th>
      <th>PARENT_ISO</th>
      <th>ISO3</th>
      <th>SUPP_INFO</th>
      <th>CONS_OBJ</th>
      <th>geometry</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>3014.0</td>
      <td>3014</td>
      <td>1</td>
      <td>Changbaishan</td>
      <td>Changbaishan</td>
      <td>UNESCO-MAB Biosphere Reserve</td>
      <td>UNESCO-MAB Biosphere Reserve</td>
      <td>International</td>
      <td>Not Applicable</td>
      <td>Not Applicable</td>
      <td>...</td>
      <td>Not Reported</td>
      <td>Not Reported</td>
      <td>Not Reported</td>
      <td>526</td>
      <td>CN-22</td>
      <td>CHN</td>
      <td>CHN</td>
      <td>Not Applicable</td>
      <td>Not Applicable</td>
      <td>POLYGON ((128.07660 42.42638, 128.09174 42.425...</td>
    </tr>
    <tr>
      <th>1</th>
      <td>3015.0</td>
      <td>3015</td>
      <td>1</td>
      <td>Wolong Nature Reserve</td>
      <td>Wolong Nature Reserve</td>
      <td>UNESCO-MAB Biosphere Reserve</td>
      <td>UNESCO-MAB Biosphere Reserve</td>
      <td>International</td>
      <td>Not Applicable</td>
      <td>Not Applicable</td>
      <td>...</td>
      <td>Not Reported</td>
      <td>Not Reported</td>
      <td>Not Reported</td>
      <td>526</td>
      <td>CN-51</td>
      <td>CHN</td>
      <td>CHN</td>
      <td>Not Applicable</td>
      <td>Not Applicable</td>
      <td>POLYGON ((103.16235 31.32249, 103.17282 31.311...</td>
    </tr>
  </tbody>
</table>
<p>2 rows × 31 columns</p>
</div>

## 3. Creating mask object, adding and manipulating layers

The mask object consists of multiple layers and manipulations performed on them. To add a layer, the four methods below perform same functions. A user may add a layer to the mask by specifying paths when a new instance is created, or use the `add_layer` method. We will add the following two files: `china_elevation.tif`, and `china_modis.tif` to the `China` mask, and name them `elevation` and `modis` layers.

```python
# Method 1: Initialize one layer, add one layer
china = geodata.Mask("China", layer_path = elevation_path)
china.rename_layer('china_elevation', 'elevation')
china.add_layer(modis_path, layer_name = 'modis')
```

```python
# Method 2:  Initialize empty, add two layers using dict
china = geodata.Mask("China")
china.add_layer(layer_path = {'elevation': elevation_path,
                             'modis': modis_path})
```

```python
# Method 3:  Initalize with two layers passed as list
china = geodata.Mask("China", layer_path = [elevation_path, modis_path],
             layer_name = ['elevation', 'modis'])
```

```python
# Method 4:  Initialize with two layers passed as dict
china = geodata.Mask("China", layer_path = {'elevation': elevation_path,
                              'modis': modis_path})
```

Display the mask object in the jupyter notebook:

```python
china
```
Mask China: 

2 layers: ['elevation', 'modis'].

No merged_mask ready. 

No shape has been extracted. 

Mask has not been saved/updated. 



Each mask object has several attributes:

- `layers`: a dictionary of name (key)  - rasterio file opener (values). The <\open DatasetReader> can be the input for many other mask methods for the module.
- `merged_mask`: the merged and flatten mask of its layers, the merged raster from `layers`
- `shape_mask`: similar to the `layers` attribute, but a dictionary of extracted shapes from the merged mask by default. Users may also extracted shape masks from specified layers in `self.layers`.
- `saved`: whether this mask object has been saved locally.
- `mask_dir`: the directory to save the mask object, by default it should be the mask dir in config.py.

Show the `slope` layer in mask `china`. The `show` method will always try to show the proper latitude and longitude, unless we call it `show(layer, lat_lon = False)`.

```python
show(china.layers['elevation'], title = 'Elevator of China in meters')
```
![png](https://github.com/east-winds/geodata/blob/mask/images/mask_creation_workflow/output_30_0.png)

**Some useful methods to examine the layers**

- `china.get_res()`: get resolution of each layer, in lat-lon coordinates
- `china.get_res(product = True)`: get grid cell size, in product of lat-lon coordinate differences
- `china.get_bounds()`: get bounds, in lat-lon coordinates

```python
china.get_bounds()
```
{'elevation': BoundingBox(left=72.99253346658085, bottom=17.994777314571174, right=136.0003674947241, top=58.23031889028454),

 'modis': BoundingBox(left=44.81844021698865, bottom=17.092253634655307, right=179.98824371508368, top=58.38850413860287)}
 
 
Note that the modis layer has a very different bounding box then the slope layer in lat-lon coordinate system. This is because the modis layer was converted to the lat-lon CRS from a different CRS when it was added to the object. The following section will explore CRS conversion.

### 3.1 CRS conversion, trimming, and cropping (if necessary)

Method `open_tif` can open a layer without adding it to the layer, this allows us to visualize it before-hand. It is a good practice to close the raster after opening it to avoid writing permission conflict issues. Closing the raster below does not involve any layer operation associated with the mask object.

```python
modis_opener = geodata.mask.open_tif(modis_path, show=True)
modis_opener.close()
```
![png](https://github.com/east-winds/geodata/blob/mask/images/mask_creation_workflow/output_36_0.png)

We can use `remove_layer` method to remove a layer to mask `china`. This method will properly close the raster file, because the raster file would remain open after being added to the mask.

```python
china.remove_layer('modis')
```
The `add_layer` method incorporates coordinate reference system (CRS) conversion to lat-lon (EPSG:4326), if necessary. Note that this method will overwrite the layer by default, if it is in the object already, unless the user specifies `replace=False`.

The method will automatically trim the all-zero columns/rows. By default, the paramater `trim` is set to `True`. If we do not set it to True, we might generate a converted raster with new CRS but many all-zero columns and rows.

```python
china.add_layer(modis_path, 'modis', trim = False)
show(china.layers['modis'], title = 'China Modis CRS converted (No trimming)')
```
![png](https://github.com/east-winds/geodata/blob/mask/images/mask_creation_workflow/output_40_1.png)

We can also crop a raster/layer with user-defined dimensions: method `crop_layer` can take either starting indices of top/left, ending indices of right/bottom, or coordinates values in lat/long to trim the raster.

The difference between `crop_layer` and `trim_layer` is that `crop_layer` must take in user specified range to crop the raster, and `trim_layer` would remove the all zero rows and columns automatically for a raster. So that if the user do not know which index to start and end to remove the empty rows/columns, `trim_raster` is better.

The method `crop_raster` (`geodata.mask.crop_raster`) is similar to `crop_layer` but can take a layer name as input, so that the user does not need to add a raster as a layer to call that method. (Similar method: `trim_layer`/`trim_raster`, `binarize_layer`/`binarize_raster`)

```python
china.crop_layer('modis', bounds = (73, 17, 135, 54))
show(china.layers['modis'], title = 'China Modis Layer Cropped')
```
![png](https://github.com/east-winds/geodata/blob/mask/images/mask_creation_workflow/output_42_0.png)

This performs the same function by passing the layer to `crop_raster`:

```python
china.layers['modis'] = geodata.mask.crop_raster(china.layers['modis'], (73, 17, 135, 54))
```
### 3.2 Filter a layer

The mask module also supports filtering a layer based on list of categorical values, a minimum (lower) boundary, or maximum (upper) boundary.

In the `filter_raster` method, a user may specify any of the `value` (the list of numberic values in the raster array to be selected), `max_bound`, and `min_bound` parameters to selected desired values. If the parameter `binarize` is False (by default), the method will return the original values of the raster that satisfy the conditions, otherwise the method will return 1 for the values that satisfy the conditions and 0 elsewhere.

#### a). Select categorical values from modis layer

Since the modis layer has 17 distinct values for different land use types, we want to create a layer of binary values, indicating unavailable land as 0, and available land as 1.

We wish to create a mask where :

- all forested areas (values 1-5) are 0 (i.e., unsuitable)
- all urban areas (13) are 0
- all others are 1

Let us use method `filter_raster` to create a layer of `modis_filtered` binary mask, where 1, 2, 3, 4, 5, and 13 will be unavailable land assigned 0 and the rest of the values will be 1 (available).

```python
avail_values = list(set(range(1, 18)) - set([1, 2, 3, 4, 5, 13]))
avail_values
```
[6, 7, 8, 9, 10, 11, 12, 14, 15, 16, 17]
```python
china.layers['modis_filtered'] = geodata.mask.filter_raster(
                china.layers['modis'], 
                binarize = True,
                values = avail_values)
```
```python
china.remove_layer('modis')
show(china.layers['modis_filtered'])
```
![png](https://github.com/east-winds/geodata/blob/mask/images/mask_creation_workflow/output_50_0.png)

#### b). Filter elevation layer

Because we cannot build renewable energy in areas with high elevation, let us set the constraint from the `elevation` layer, by using elevation < 4000m at 1 and other areas as 0. The result layer `elevation_filtered` will have only 1 and 0 as unique values.

```python
china.filter_layer('elevation', 
                   dest_layer_name = 'elevation_filtered', 
                   max_bound = 4000,
                   binarize = True)
```
```python
china.remove_layer('elevation')
show(china.layers['elevation_filtered'])
```
![png](https://github.com/east-winds/geodata/blob/mask/images/mask_creation_workflow/output_53_0.png)

#### c). Filter slope layer

We also cannot build renewable energy in area with large slopes, so let us set another constraint from the `slope` layer from the slope tif file, by using slope < 20 degree at 1 and else as 0. The result layer `slope_filtered` will have only 1 and 0 as unique values.

</div>

First, add the slope raster to the china mask.

```python
china.add_layer(slope_path, layer_name = 'slope')
show(china.layers['slope'])
```
![png](https://github.com/east-winds/geodata/blob/mask/images/mask_creation_workflow/output_56_1.png)

Filter the raster, delete the old slope layer.

```python
china.filter_layer('slope', 
                   dest_layer_name = 'slope_filtered', 
                   max_bound = 20,
                   binarize = True)
```
```python
china.remove_layer('slope')
show(china.layers['slope_filtered'])
```
![png](https://github.com/east-winds/geodata/blob/mask/images/mask_creation_workflow/output_59_0.png)

### 3.3 Additional show options

We can plot the provinces on a selected layer by taking `shape` input in the `show()` method. Here, we will use the `china_shapes` that we obtained from `all_shape`. Its `geometry` column is a Series of shapes (shapely.geometry or MultiPolygon) for Chinese provinces.

```python
show(china.layers['modis_filtered'], shape = china_shapes['geometry'])
```
![png](https://github.com/east-winds/geodata/blob/mask/images/mask_creation_workflow/output_62_0.png)

## 4. Opening a shapefile and adding shape features as a layer

Recall that we have previously loaded the environmental protected shapes of China in a GeoPandas dataframe.

```python
len(wdpa_shapes)
```
78



The three shapefiles have 78 features altogether, but we want to add all the features to one new layer instead of 78 new layers. The input shape should be a python dictionary, where there is a key for each unique shape. Also, in the `add_shape_layer` method, we will specify a `combine_name` to combine the features into one layer in this case, since we want the mask to have just one more layers, not 78 more layers.

When adding a shapefile, we must specify the dimensions. We will also use `reference layer = 'slope_filtered'` so the new shape layer will have the same dimension with the `slope_filtered` layer. If the mask is empty and does not contain any layer, the user will have to specify the `resolution` parameter for the raster layer dimension.

By default, this method will have paramater `exclude` that defaults to `False`. When it is true, area inside the shape is 0. When it is false, area inside the shape is 1. In this use case, however, we want 0 for area inside of the shape as they are environmental protected areas to exclude. We can just use the default method call.

```python
china.add_shape_layer(wdpa_shapes['geometry'].to_dict(), 
                      reference_layer = 'slope_filtered', 
                      combine_name = 'protected')
show(china.layers['protected'], 
     title = 'WDPA Protected area shape features as a new layer',
     grid = True)
```
![png](https://github.com/east-winds/geodata/blob/mask/images/mask_creation_workflow/output_67_1.png)

We can also use the parameter `buffer` in `add_shape_layer` method to create an approximate representation of all locations within a given (perpindicular) distance of the shape object. The units for the buffer are given in kilometers.

Note that since the units of the original shape are in lat-lon coordinates, when we add the buffer, we will need to have a CRS that has meter as unit. The program will convert the shapes to that CRS, add the buffer around shapes, then convert it back to the lat-lon CRS system. By default, we used "EPSG:6933", an equal area projection CRS to add buffer in kilometer.

```python
km_buffer = 20

china.add_shape_layer(wdpa_shapes['geometry'].to_dict(), 
                      reference_layer = 'slope_filtered', 
                      combine_name = 'protected_with_buffer',
                      buffer = km_buffer)

show(china.layers['protected_with_buffer'], 
     title = f"WDPA Protected area shape with {km_buffer}km buffer",
     grid = True)

china.remove_layer('protected_with_buffer')
```
![png](https://github.com/east-winds/geodata/blob/mask/images/mask_creation_workflow/output_69_1.png)

## 5. Merging and flattening layers

In order to combine all layers into one, we use the `merge_layer` method which creates a new layer called `merged_mask`. This merges multiple layers together and flattens them using either **and** (default) or **sum** method, saving the result to `self.merged_mask` by default. Geospatial bounds and resolution of the output layer are in the units of the input file coordinate reference system, but by default, we will use the resolution of the layer with the best (finest) resolution for the output bounds/resolution, unless a reference layer is provided. In this case, the resolution of the merged_mask is the same with the `modis_filtered` layer.

```python
china.get_res()
```
{'modis_filtered': (0.006363926718366056, 0.006364039220827179),
 'elevation_filtered': (0.008983152841195215, 0.008983152841195215),
 'slope_filtered': (0.008983152841195215, 0.008983152841195215),
 'protected': (0.008983152841195215, 0.008983152841195215)}
```python
china.merge_layer(attribute_save = False, show = False).res
```
(0.006363926718366056, 0.006364039220827179)
### 5.1 binary AND method

By default, the `merge_layer` method will use a binary 'and' method: for each grid cell, if any of the n layers are 0, then the returned `self.merged_layer` will also have 0 at that location. In other words, if all the layers indicate that a land is available (!=0), the merged result will have value 1.

`merge_layer` may also take in an optional parameter `layers`, which is a list of layer names stored in the object, if the user does not wish to merge all layers in the object. If the user does not want to save the result to the `merged_mask` attribute, the user can specify `attribute_save = False`.

```python
# merge and plot only, do not save
china.merge_layer(attribute_save = False, layers = ['slope_filtered', 'modis_filtered'])
```
![png](https://github.com/east-winds/geodata/blob/mask/images/mask_creation_workflow/output_76_0.png)

<open DatasetReader name='/vsimem/47a585bd-eec9-4672-8ee7-2e20d120ccda/47a585bd-eec9-4672-8ee7-2e20d120ccda.tif' mode='r'>
  
```python
china
```
  
Mask China: 
  
4 layers: ['modis_filtered', 'elevation_filtered', 'slope_filtered', 'protected'].
  
No merged_mask ready. 
  
No shape has been extracted. 
  
Mask has not been saved/updated. 
  
Try again with the `reference_layer` parameter:
  

```python
china.merge_layer(layers = ['elevation_filtered', 'modis_filtered'], reference_layer = 'elevation_filtered', show = False)
```

The result of the `merged_mask` method is saved to `china.merged_mask` with the same resolution as the reference layer, in this case `elevation_filtered`.

```python
china.merged_mask.res
```
(0.008983152841195215, 0.008983152841195215)
  
  
For the purpose of this demonstration, we will select the AND method for the final merged_mask. We can also trim the border of the merged mask since the 4 layers have different boundaries. We can set the parameter `trim = True`.

```python
china.merge_layer(trim = True)
```
![png](https://github.com/east-winds/geodata/blob/mask/images/mask_creation_workflow/output_83_0.png)


### 5.2 sum method

The sum method will add up the values from all the layers using weights. When there is no weight dict provided, all the layers for merging will have weights of 1 by default.

Note: since we are not using the sum method to proceed to the following sections, we will keep `attribute_save = False` to prevent this method from overwriting the mask we have previously created above.

```python
china.merge_layer(method = 'sum', 
                  attribute_save = False,
                  trim = True)
```

![png](https://github.com/east-winds/geodata/blob/mask/images/mask_creation_workflow/output_86_1.png)
  

<open DatasetReader name='/vsimem/1bab6b72-616b-4f7a-995a-9edb7c795457/1bab6b72-616b-4f7a-995a-9edb7c795457.tif' mode='r'>
  
This distribution is completely arbitrary for the purpose of demonstration of the module: (Note: The weights do not need to have a total of 1)

- elevation_filtered: 0.15, slope_filtered: 0.1, modis_filtered: 0.3, protected: 0.45

We will write the result to a new variable `customized_merged_layer` for continuing processing.

```python
customized_merged_layer = china.merge_layer(method = 'sum', weights = {
        'elevation_filtered': 0.15,
        'slope_filtered': 0.1,
        'modis_filtered': 0.3,
        'protected': 0.45
    }, attribute_save = False, trim = True)
```
![png](https://github.com/east-winds/geodata/blob/mask/images/mask_creation_workflow/output_88_0.png)

If the continuous value created by `merged_mask` represents a suitability metric, we could set a minimum value of 0.8 to be considered "suitable" (or 1). We then apply the `filter_raster` method on the merged layer.

```python
customized_merged_layer = geodata.mask.filter_raster(customized_merged_layer, min_bound = 0.8, binarize = True)
show(customized_merged_layer)
```
![png](https://github.com/east-winds/geodata/blob/mask/images/mask_creation_workflow/output_90_0.png)

## 6. Eliminate small contiguous areas

Using the above methods, we might end up with many small contiguous areas that are marked suitable but surrounded by an unsuitable region. We may want to exclude such regions from renewable energy development. The `filter_area` method will remove the small contiguous suitable regions by transforming the merged mask raster to polygons/shapes, calculating the area of each polygon, and filtering out polygons that are smaller than a given threshold. Units are given in kilometer-squared (km$^2$).

By default, `filter_area` uses the merged mask raster and returns a new raster, unless input/output layers are specified by `layer_name` and `dest_layer_name`.

By default, its `shape_value` parameter is 1, indicating that we are only interested in finding all groups of cells with value 1 (suitable) for elimination. We specify the threshold with the `min_area` parameter.

Note: the `filter_area` method may take a long time (5 or more minutes depending on the complexity of your layer and your computational setup). The method relies upon `rasterio.rasterize`, see performance notes: https://rasterio.readthedocs.io/en/latest/api/rasterio.features.html#rasterio.features.rasterize

For example, if we focus on Guangdong province in Southern China from the merged mask, we notice that there are many small islands in the ocean that are marked as suitable areas. We want to exclude these small regions from our merged mask.

```python
plt.imshow(china.merged_mask.read(1)[4800:5300, 5700:6600], interpolation = 'none')
plt.show()
```
![png](https://github.com/east-winds/geodata/blob/mask/images/mask_creation_workflow/output_94_0.png)

Call `filter_area` to remove all contiguous suitable region shapes smaller than 100 km$^2$:

```python
china.merged_mask = geodata.mask.filter_area(china, min_area = 100)
```

There shapes are removed in the new merged_mask.

```python
plt.imshow(china.merged_mask.read(1)[4800:5300, 5700:6600], interpolation = 'none')
plt.show()
```
![png](https://github.com/east-winds/geodata/blob/mask/images/mask_creation_workflow/output_98_0.png)

## 7. Extracting shapes from mask

Sometimes the user needs to generate masks and perform analysis for a collection of regions (e.g., at the state/province level). The purpose of shape extraction (`extract_shapes`) is to separate `merged_mask` values for each region, with the result a dictionary of name-mask pairs in the `shape_mask` attribute of the mask object. The values of `shape_mask` will be 0 outside of the shape, and will be `merged_mask` inside of the shape.

For the purpose of this demonstration, we will only select the province of Jiangsu, Zhejiang, and Shanghai.

```python
china_shapes_subset = china_shapes[china_shapes['name'].isin(
                                    ['Jiangsu', 'Zhejiang', 'Shanghai'])]
china_shapes_subset
```
<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>featurecla</th>
      <th>scalerank</th>
      <th>adm1_code</th>
      <th>diss_me</th>
      <th>iso_3166_2</th>
      <th>wikipedia</th>
      <th>iso_a2</th>
      <th>adm0_sr</th>
      <th>name</th>
      <th>name_alt</th>
      <th>...</th>
      <th>name_nl</th>
      <th>name_pl</th>
      <th>name_pt</th>
      <th>name_ru</th>
      <th>name_sv</th>
      <th>name_tr</th>
      <th>name_vi</th>
      <th>name_zh</th>
      <th>ne_id</th>
      <th>geometry</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>2113</th>
      <td>Admin-1 scale rank</td>
      <td>2</td>
      <td>CHN-1820</td>
      <td>1820</td>
      <td>CN-33</td>
      <td>None</td>
      <td>CN</td>
      <td>5</td>
      <td>Zhejiang</td>
      <td>Zhèjiāng</td>
      <td>...</td>
      <td>Zhejiang</td>
      <td>Zhejiang</td>
      <td>Zhejiang</td>
      <td>Чжэцзян</td>
      <td>Zhejiang</td>
      <td>Zhejiang</td>
      <td>Chiết Giang</td>
      <td>浙江省</td>
      <td>1159311795</td>
      <td>MULTIPOLYGON (((121.25553 28.14000, 121.26454 ...</td>
    </tr>
    <tr>
      <th>2114</th>
      <td>Admin-1 scale rank</td>
      <td>2</td>
      <td>CHN-1819</td>
      <td>1819</td>
      <td>CN-31</td>
      <td>None</td>
      <td>CN</td>
      <td>5</td>
      <td>Shanghai</td>
      <td>Shànghǎi</td>
      <td>...</td>
      <td>Shanghai</td>
      <td>Szanghaj</td>
      <td>Xangai</td>
      <td>Шанхай</td>
      <td>Shanghai</td>
      <td>Şanghay</td>
      <td>Thượng Hải</td>
      <td>上海市</td>
      <td>1159311855</td>
      <td>MULTIPOLYGON (((121.83774 31.37515, 121.85377 ...</td>
    </tr>
    <tr>
      <th>2116</th>
      <td>Admin-1 scale rank</td>
      <td>2</td>
      <td>CHN-1818</td>
      <td>1818</td>
      <td>CN-32</td>
      <td>None</td>
      <td>CN</td>
      <td>3</td>
      <td>Jiangsu</td>
      <td>Jiāngsū</td>
      <td>...</td>
      <td>Jiangsu</td>
      <td>Jiangsu</td>
      <td>Jiangsu</td>
      <td>Цзянсу</td>
      <td>Jiangsu</td>
      <td>Jiangsu</td>
      <td>Giang Tô</td>
      <td>江苏省</td>
      <td>1159311859</td>
      <td>MULTIPOLYGON (((119.88071 32.12519, 119.90138 ...</td>
    </tr>
  </tbody>
</table>
<p>3 rows × 84 columns</p>
</div>

Converting the filtered shape dictionary to a python dictionary as the input for `extract_shapes`, where the keys for the dictionary will be the names of the new extracted shape layers.

```python
china_shapes_subset = china_shapes_subset[['name', 
                                           'geometry']].set_index("name")['geometry'].to_dict()
china_shapes_subset
```
{'Zhejiang': <shapely.geometry.multipolygon.MultiPolygon at 0x26f6eb5fbe0>,
  
 'Shanghai': <shapely.geometry.multipolygon.MultiPolygon at 0x26f6eb5fbb0>,
  
 'Jiangsu': <shapely.geometry.multipolygon.MultiPolygon at 0x26f6eb5fca0>}
  
Extract the shapes from the merged_mask.

```python
china.extract_shapes(china_shapes_subset)
```
INFO:geodata.mask:Extracted shape Zhejiang added to attribute 'shape_mask'.
INFO:geodata.mask:Extracted shape Shanghai added to attribute 'shape_mask'.
INFO:geodata.mask:Extracted shape Jiangsu added to attribute 'shape_mask'.
The resulting mask object contains the dictionary `shape_mask` with the extracted values:

```python
china
```
Mask China: 
  
4 layers: ['modis_filtered', 'elevation_filtered', 'slope_filtered', 'protected'].
  
Merged_mask merged/flattened. 
  
3 shape_mask: ['Zhejiang', 'Shanghai', 'Jiangsu']. 
  
Mask has not been saved/updated. 
  
## 8. Saving and Loading masks

```python
china.save_mask()
```
With the mask saved, the user can now load the layers or shapes with `xarray` instead if preferred.

```python
shape_xr_lst = china.load_shape_xr()
shape_xr_lst['Zhejiang'].plot()
```





<matplotlib.collections.QuadMesh at 0x26f11e347c0>
![png](https://github.com/east-winds/geodata/blob/mask/images/mask_creation_workflow/output_112_2.png)

Optional: closing all the files when saving the mask. This can avoid possible write permission error.

```python
china.save_mask(close_files = True)
```
  
Loading a previously saved mask.

```python
china_2 = geodata.mask.load_mask("china")
```

```python
china_2
```
Mask china: 
  
4 layers: ['elevation_filtered', 'modis_filtered', 'protected', 'slope_filtered'].
  
Merged_mask merged/flattened. 
  
3 shape_mask: ['Jiangsu', 'Shanghai', 'Zhejiang']. 
  
Mask has been saved. 
  

