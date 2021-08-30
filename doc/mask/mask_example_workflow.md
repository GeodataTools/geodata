# Mask Workflow Example

A short guide to use the mask module to process geospatial files and extract on cutouts.

For a complete jupyter notebook example, see: [mask.ipynb](https://github.com/east-winds/geodata/tree/mask/tests/mask_test.ipynb)

Key methods described here (WITH LINKS):

- [add_layer](#3.-Creating-mask-object,-adding-and-manipulating-layers)

## 1. Introduction

Geodata is able to process geospatial data to extract cutouts over specified geographies. Built off the [rasterio library](https://rasterio.readthedocs.io/en/latest/quickstart.html), the mask module imports rasters and shapefiles, merges and flattens multiple layers together, and extracts subsetted cutout data from merged masks and shapefiles.

Its current functionalities in details are:

- Creating mask object, adding and manipulating layers
- Openning shapefile and adding shape features as a layer
- Merging and flattening layers
- Extracting shapes from mask
- Saving and Loading masks

## 2. Setup

To start, import the required dependencies:

```
import geodata
import numpy as np
from geodata.mask import show
```

To launch a logger for detailed debugging, run:

```
import logging
logging.basicConfig(level=logging.INFO)
```

We use [cartopy](https://scitools.org.uk/cartopy/docs/latest/tutorials/using_the_shapereader.html#cartopy.io.shapereader.Reader) to download some common administrative region shapes, but user-provided shapefiles will also work:

```
import cartopy.io.shapereader as shpreader
```

## 3. Creating mask object, adding and manipulating layers

The mask object consists of multiple layers and manipulations performed on them. To add a layer, the four methods below perform same functions. A user may add a layer to the mask by specifying paths when a new instance is created, or use the `add_layer` method. We will add the following two files: `elevation_slope.tif`, and `MODIS_China.tif` to the `China` mask, and name them `slope` and `modis` layers.

```
china = geodata.Mask("China", layer_path = slope_path)
china.rename_layer('elevation_slope', 'slope')
china.add_layer(modis_path, layer_name = 'modis')
```

```
china = geodata.Mask("China")
china.add_layer(layer_path = {'slope': slope_path,
                             'modis': modis_path})
```

```
china = geodata.Mask("China", layer_path = [slope_path, modis_path],
             layer_name = ['slope', 'modis'])
```

```
china = geodata.Mask("China", layer_path = {'slope': slope_path,
                              'modis': modis_path})
```

Display the mask object in the jupyter notebook:

```
>>> china
Mask China: 
2 layers: ['slope', 'modis'].
No merged_mask ready. 
No shape has been extracted. 
Mask has not been saved/updated. 
```

Each mask object has several attributes:

- `layers`: a dictionary of name (key)  - rasterio file opener (values). The <\open DatasetReader> can be the input for many other mask methods for the module.
- `merged_mask`: the merged and flatten mask of its layers, the merged raster from `layers`
- `shape_mask`: similar to the `layers` attribute, but a dictionary of extracted shapes from the merged mask by default. Users may also extracted shape masks from specified layers in `self.layers`.
- `saved`: whether this mask object has been saved locally.
- `mask_dir`: the directory to save the mask object, by default it should be the mask dir in config.py.

Show the `slope` layer in mask `china`. The `show` method will always try to show the proper latitude and longitude.

```
show(china.layers['slope'])
```

![slope](https://github.com/east-winds/geodata/tree/mask/images/mask_test/slope.png)

Some useful methods for the layers:

- `china.get_res()`: get resolution of each layer, in lat-lon coordinates
- `china.get_res(product = True)`: get grid cell size, in product of lat-lon coordinate differences
- `china.get_bounds()`: get bounds, in lat-lon coordinates

```
>>> china.get_bounds()
{'slope': BoundingBox(left=72.99253346658085, bottom=17.994777314571174, right=136.0003674947241, top=58.23031889028454),
 'modis': BoundingBox(left=44.81844021698865, bottom=17.092253634655307, right=179.98824371508368, top=58.38850413860287)}
```

Note that the modis layer has a very different bounding box then the slope layer in lat-lon coordinate system. This is because the modis layer was converted to the lat-lon CRS from a different CRS when it was added to the object. The following section will explore the behind scene.

### 3.1 CRS conversion, trimming, cropping (if necessary)

Method `open_tif` can open a layer without adding it to the layer, this allows us to visualize it before-hand. It is a good practice to close the raster after openning it to avoid writing permission conflict issues.

```
modis_openner = geodata.mask.open_tif(modis_path, show=True)
modis_openner.close()
```

The `add_layer` method incorporates CRS conversion, let us see what the layer looks like after we add it to the object.

We can use `remove_layer` method to remove a layer to mask `china`. This method will properly close the raster file.

```
china.remove_layer('modis')
```

The `add_layer` method incorporates CRS conversion. When the user pass in a path to a raster file, the program will automatically convert the CRS of the raster to the lat-lon CRS by default.

Note that this method will overwrite the layer by default, if it is in the object already, unless the user specifies `replace=True`.

The method will automatically trim the all-zero columns/rows. By default, the paramater `trim` is set to `True`. If we do not set it to True, we might generate a converted raster with new CRS but many all-zero columns and rows.

```
china.add_layer(modis_path, 'modis', trim = False)
show(china.layers['modis'], title = 'China Modis CRS converted (No trimming)')
```

We can also **arbitrary** crop a raster/layer: method `crop_layer` can take either starting indices of top/left, ending indices of right/bottom, or coordinates values in lat/long to trim the raster.

We also have a method `crop_raster` (`geodata.mask.crop_raster`) similarly to `crop_layer` but we can have any raster as input, which indicates that users do not need to add a raster as a layer to call that method. (Similar method: `trim_layer`/`trim_raster`, `binarize_layer`/`binarize_raster`)

```
china.crop_layer('modis', bounds = (73, 17, 135, 54))
show(china.layers['modis'], title = 'China Modis Layer Cropped')
```

The line below performs the same function with the code block above in arbitrary cropping the modis layer.

```
china.layers['modis'] = geodata.mask.crop_raster(china.layers['modis'], (73, 17, 135, 54))
```

### 3.2 Categorical Value Extraction, if necessary

However, the modis layer have 17 distinct values, we may want to create a layer of binary values, indicating unavailable land as 0, and available land as 1.

Values 1, 2, 3, 4, 5 are 5 types of forest for the modis layer, let us use method `binarize_raster` to create a layer of `modis_forest` mask with 1 and 0. 1, 2, 3, 4, 5 will be unavailable land, therefore we need to take in the rest of the values to make them 1 (available).

```
values = np.arange(6, 18)
china.layers['modis_forest'] = geodata.mask.binarize_raster(china.layers['modis'], values = values)
china.remove_layer('modis')
show(china.layers['modis_forest'])
```

## 4. Openning shapefile and adding shape features as a layer

Check attributes in the shapes contained in path `enep_shape_path`. The `shape_attribute` method checks the attribute of the features in the shapefile by showing the attribute-value pair of the first feature in the shapefile.

```
geodata.mask.shape_attribute(enep_shape_path)
```

Loading shapes with the `get_shape` method. It is more convenient to load the shapes into a python dictionary as a input for `add_shape_layer` method. If we do not specify `return_dict = True`, the result will be a geopandas dataframe. We will add the shapes to the china mask object.

There are more examples of the `get_shape` method later in this demo when we want to extract the chinese provinces from the mask.

```
>>> protected_area_shapes = geodata.mask.get_shape(enep_shape_path, key = 'NAME', return_dict = True)
>>> len(protected_area_shapes)
824
```

There are 824 shapes, but we want to add all the shapes to one new layer instead of 824 new layers. Therefore, in the `add_shape_layer` method, we will pecify `combine_name` so that the program will combine the features to one layer with the `combine_name` as its layer name.

We will also use `reference layer = 'slope'` so the new shape layer will have the same dimension with the `slope` layer. If the mask is empty and does not contain any layer, the user will have to specify the dimension as another parameter.

By default, this method will have paramater `invert` that defaults to `True`, and generate 1 for area covered by the shape, and 0 for the area outside of the shape. In this use case, however, we want 0 for area inside of the shape as they are the environmental protected area. We will specify `invert = False`.

```
china.add_shape_layer(protected_area_shapes, 
                      reference_layer = 'slope', 
                      combine_name = 'protected',
                      invert = False)
show(china.layers['protected'], title = 'Protected area shape features as a new layer')
```

## 5. Merging and flattening layers

`merge_layer` method is important for merging all the layers, or the specified layers, into one single mask in object attribute called `merged_mask`.

It merges multiple layers together and flatten it using either **and** (default) or **sum** method, saving the result to `self.merged_mask` by default. Geospatial bounds and resolution of a new output file are in the units of the input file coordinate reference system, but by default, we will use the resolution of the layer with the best (finest) resolution for the output bounds/resolution, unless a reference layer is provided. In this case, the resolution of the merged_mask is the same with the `modis_forest` layer.

```
>>> china.get_res()
{'slope': (0.008983152841195215, 0.008983152841195215),
 'modis_forest': (0.006363926718366056, 0.006364039220827179),
 'protected': (0.008983152841195215, 0.008983152841195215)}
```

```
>>> china.merge_layer(attribute_save = False, show = False).res
(0.008983152841195215, 0.008983152841195215)
```

### 5.1 binary AND method

By default, the `merge_layer` method will use a binary 'and' method: if any of the n grid cells of the n layers at the same location have 0, then the returned `self.merged_layer` will also have 0 at that location. In other words, if all the layers indicate that a land is not unavailable (!=0), the merged result will have value 1.

`merge_layer` may also take in an optional parameter `layers`, which is a list of layer names stored in the object. If users do not wish to create the final merged mask with all layers, they can specify which layers to use. If the user do not want to save the result of this method to the `merged_mask` attribute, the user can specify `attribute_save = False`.

```
china.merge_layer(attribute_save = False, layers = ['slope', 'modis_forest'])
```

Try a different reference layer.

```
china.merge_layer(layers = ['slope', 'modis_forest'], reference_layer = 'slope', show = False)
```

Now the result of the `merged_mask` method is saved to `china.merged_mask`, and the resolution of the `merged_mask` is the same with the `slope` reference layer.

```
>>> china.merged_mask.res
(0.008983152841195215, 0.008983152841195215)
```

### 5.2 sum method

The sum method will add up the values from all the layers. We can also customize the weights. The behind scene of this method is that it multiplys each layers with the corresponding weight, and add the in-memory temporary layers together.

```
china.merge_layer(method = 'sum')
```

This distribution is completely arbitrary for the purpose of demonstration of the module:

- slope: 25%, modis_forest: 30%, protected 45%

The weights do not need to have a total of 1.

```
china.merge_layer(method = 'sum', weights = {
    'slope': 0.25,
    'modis_forest': 0.3,
    'protected': 0.45
})
```

We can also trim the border of the merged mask since not 4 layers have the same boundary, and the border values are not useful. We can set the parameter `trim = True`.

```
china.merge_layer(method = 'sum', weights = {'slope': 0.25, 'modis_forest': 0.3, 'protected': 0.45},
                 trim = True)
```

## 6.0 Extracting shapes from mask

Extracting shapes from mask is different from section 4, where we added shape features as a new layer. The purpose of shape extraction is that the user may need to generate masks and perform analysis on a state/province level. After generating the `merged_mask` with appropriate mask values, the program can generate the `merged_mask` for every state or province. The result of shape extraction is a dictionary of state_name - state mask pair in `shape_mask` attribute of the mask object. It is still the case that the values in the mask will be 0 outside of the shape in the mask, but inside of the shape of the province/state, we will have the merged_mask values.

Let us get province shapes from `cartopy` and save the path as `prov_path`. This can also be the path to user-supplied shape files.

```
prov_path = shpreader.natural_earth(resolution='10m', category='cultural', name = 'admin_1_states_provinces')
prov_path
```

Check attributes in the shapes contained in path `prov_path`

```
geodata.mask.shape_attribute(prov_path)
```

For the `get_shape` method, while the targets list the exact names/value of key_name attributes of the shapes, users may also ignore it and use condition_key and condition_value to find the desired shapes. For example, the call below will find all the shapes of provinces that belongs to China:
`mask.get_shape(path_to_province_shapes, key = 'name', condition_key = 'admin', condition_value = 'China')`

```
china_all_shapes = geodata.mask.get_shape(prov_path, key = 'name_en',
condition_key = 'admin', condition_value = 'China')
china_all_shapes
```

We can also ignore condition, just take three provinces of China by naming them out:

```
china_shapes = geodata.mask.get_shape(prov_path, key = 'name_en',
                         targets = ['Jiangsu', 'Zhejiang', 'Shanghai'],
                         return_dict = True)
```

Extract the shapes from the merged_mask.

```
>>> china.extract_shapes(china_shapes)
>>> china
Mask China:
4 layers: ['bins', 'forest', 'slope', 'modis_forest'] .
Merged_mask merged/flattened.
3 shape_mask: ['Zhejiang', 'Shanghai', 'Jiangsu'].
Mask has not been saved/updated.
```

## 7. Saving and Loading masks

```
>>> china.save_mask()
INFO:geodata.mask:Mask China successfully saved at D:/Users/davison_lab_data/masks
```

Note that since "Mask has been saved", we can now load the layers or shapes with xarray.

```
shape_xr_lst = china.load_shape_xr()
shape_xr_lst['Zhejiang'].plot()
```

Optional: closing all the files when saving the mask. This can avoid possible write permission error.

```
china.save_mask(close_files = True)
```

Loading a previously saved mask

```
>>> china_2 = geodata.mask.load_mask("china")
>>> china_2
Mask china:
7 layers: ['bins', 'forest', 'Jiangsu', 'modis_forest', 'Shanghai', 'slope', 'Zhejiang'] .
Merged_mask merged/flattened.
3 shape_mask: ['Jiangsu', 'Shanghai', 'Zhejiang'].
Mask has been saved.
```
