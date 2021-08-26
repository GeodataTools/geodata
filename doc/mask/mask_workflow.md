# Mask Workflow Example

A short guide to use the mask module to process geospatial files and extract on cutouts.

For a complete jupyter notebook example, see: [mask.ipynb](https://github.com/east-winds/geodata/tree/mask/tests/mask_test.ipynb)

Key methods described here (WITH LINKS):
- ...



## Introduction

Geodata is able to process geospatial data to extract cutouts over specified geographies. Built off the [rasterio library](https://rasterio.readthedocs.io/en/latest/quickstart.html), the mask module imports rasters and shapefiles, merges and flattens multiple layers together, and extracts subsetted cutout data from merged masks and shapefiles.

Its current functionalities in details are:

- Creating mask, adding layers from .tif files
- CRS conversion, cropping, trimming, binarizing layers
- Merging and flattening layers
- Adding shape files as layers
- Extracting shapes from mask layers

## Setup

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


## Creating mask object, adding and manipulating layers

The mask object consists of multiple layers and manipulations performed on them. To add a layer...

The four methods below perform same functions. A user may add a layer to the mask by specifying paths when a new instance is created, or use the `add_layer` method.


(LET'S ADD THESE LAYERS INSTEAD:
- Elevation_Slope.tif
- MODIS_China.tif
)

We will add the following two files: `FINAL_GRID_5BINS.tif`, and `FINAL_GRID_FOREST_MED.tif` to the `China` mask, and name them `bins` and `forest` layers.



```
china = geodata.mask("China", layer_path = grid_bins_path)
china.rename_layer('FINAL_GRID_5BINS', 'bins')
china.add_layer(forest_med_path, layer_name = 'forest')
```

```
china = geodata.mask("China")
china.add_layer(layer_path = {'bins': grid_bins_path,
                             'forest': forest_med_path})
```

```
china = geodata.mask("China", layer_path = [grid_bins_path, forest_med_path],
             layer_name = ['bins', 'forest'])
```

```
china = geodata.mask("China", layer_path = {'bins': grid_bins_path,
                              'forest': forest_med_path})
```

Display the object in the jupyter notebook:

```
>>> china
Mask China:
2 layers: ['bins', 'forest'].
No merged_mask ready.
No shape has been extracted.
Mask has not been saved/updated.
```

Each mask object has several attributes:
- `layers`: a dictionary of name (key)  - rasterio file opener (values). The <\open DatasetReader> can be the input for many other mask methods for the module.  


(CAN WE DELETE THIS OR MOVE THIS TO ANOTHER SECTION ON FILE I/O?)
[FYI: not important for its normal usage, but it may be helpful]: Although we save openning layerd in Rasterio, which only reads/opens/writes data on disk, sometimes we may have in-memory tif files in the `layers` attributes, and the product `merged_mask` and `shape_mask` are also in-memory temporary files. This is because we might create a new layer after automatic CRS conversion if we detect that the input file is not in latitude-longitude CRS, and in rasterio, a CRS conversion, merging-flattening, shapes on raster, cropping, and many other methods that make changes to the raster will requires creating a new file on disk, but we can avoid creating too many temporarily files and deleting them later by using memory files. Read more here: https://rasterio.readthedocs.io/en/latest/topics/memory-files.html

- `merged_mask`: the merged and flatten mask of its layers, the merged raster from `layers`
- `shape_mask`: similar to the `layers` attribute, but a dictionary of extracted shapes from the merged mask by default. Users may also extracted shape masks from specified layers in `self.layers`.
- `saved`: whether this mask object has been saved locally.
- `mask_dir`: the directory to save the mask object, by default it should be the path in config.py.

To finish adding all the layers: add another layers `Elevation_Slope.tif` as `slope`.

```
china.add_layer(slope_path, layer_name = 'slope')
```

Some useful methods for the layers:
- `china.get_res()`: get resolution of each layer, IN WHAT UNITS...
- `china.get_res(product = True)`: get grid cell size, IN WHAT UNITS
- `china.get_bounds()`: get bounds, IN WHAT UNITS...


### CRS conversion, trimming, cropping (if necessary)

Method `open_tif` can open a layer without adding it to the layer, this allows us to visualize it before-hand.

```
china.open_tif(modis_china_path, show=True)
```

The `add_layer` method incorporates CRS conversion, let us see what the layer looks like after we add it to the object.

This method will overwrite the layer by default, if it is in the object already.

```
china.add_layer(modis_china_path, layer_name = 'modis')
```

Note that now the `show` method will show the proper latitude and longitude now, once the layer is in correct CRS. However, a problem with the result below is that there are too many surrounding columns, or maybe rows, that are all zero, and we might want to trim it as a result. We can add the layer again

We can use `remove_layer` method to remove a layer to mask `china`, method `add_layer` by default can simply replace the old layer as well.

```
#same with: geodata.mask.show(china.layers['modis'])
show(china.layers['modis'])
```

(LET'S FOCUS ON THIS, THE MOST COMMON USE CASE)
If you plot this new layer, you would realize that there is too many empty columns as border, and we need to trim it while keeping all the valid data points.

```
#china.remove_layer('modis')
china.add_layer('MODIS_China.tif', layer_name = 'modis', trim_raster = True)
```

```
#china.remove_layer('modis')
china.add_layer('MODIS_China.tif', layer_name = 'modis', trim_raster = True)
```

We can also **arbitrary** crop a raster/layer: method `crop_layer` can take either starting indices of top/left, ending indices of right/bottom, or coordinates values in lat/long to trim the raster.

We also have a method `crop_raster` (`geodata.mask.crop_raster`) similarly to `crop_layer` but we can have any raster as input, which indicates that users do not need to add a raster as a layer to call that method. (Similar method: `trim_layer`/`trim_raster`, `binarize_layer`/`binarize_raster`)

```
china.crop_layer('modis', bounds = (73, 17, 135, 54))
#same thing with: china.layers['modis'] = geodata.mask.crop_raster(china.layers['modis'], (73, 17, 135, 54)) #similar to 5bins
show(china.layers['modis'])
```

### Categorical Value Extraction, if necessary

However, the modis layer have 17 distinct values, we may want to create a layer of binary values, indicating unavailable land as 0, and available land as 1.

Values 1, 2, 3, 4, 5 are 5 types of forest for the modis layer, let us use method `binarize_raster` to create a layer of `modis_forest` mask with 1 and 0. 1, 2, 3, 4, 5 will be unavailable land, therefore we need to take in the rest of the values to make them 1 (available).


```
values = np.arange(6, 18)
china.layers['modis_forest'] = geodata.mask.binarize_raster(china.layers['modis'], values = values)
china.remove_layer('modis')
```

## Merging and flattening layers

`merge_layer` method is important for merging all the layers, or the specified layers, into one single mask.

It merges multiple layers together and flatten it using either **and** (default) or **sum** method, saving the result to `self.merged_mask` by default. Geospatial bounds and resolution of a new output file are in the units of the input file coordinate reference system, but by default, we will use the resolution of the layer with the best (finest) resolution for the output bounds/resolution, unless a reference layer is provided.

#### binary AND method

By default, the `merge_layer` method will use a binary 'and' method: if any of the n grid cells of the n layers at the same location have 0, then the returned `self.merged_layer` will also have 0 at that location. In other words, if all the layers indicate that a land is not unavailable (!=0), the merged result will have value 1.

`merge_layer` may also take in an optional parameter `layers`, which is a list of layer names stored in the object. If users do not wish to create the final merged mask with all layers, they can specify which layers to use.

```
china.merge_layer(attribute_save = False, layers = ['bins', 'forest'])
```

#### sum method
(MAYBE INCLUDE A REFERENCE LAYER EXAMPLE)

The sum method will add up the values from all the layers. We can also customize the weights. The behind scene of this method is that it multiplys each layers with the corresponding weight, and add the in-memory temporary layers together.

```
china.merge_layer(method = 'sum')
```

This distribution is completely arbitrary for the purpose of demonstration of the module:
- bins: 5%, forest: 25%, slope 40%, and modis_forest 30%.

```
china.merge_layer(method = 'sum', weights = {
    'bins': 0.05,
    'forest': 0.25,
    'slope': 0.4,
    'modis_forest': 0.3
})
```

The result mask use the grid cell resolution with the layer with the **finest resolution**. We can also trim the border of the merged mask since not 4 layers have the same boundary, and the border values are not useful. We can set the parameter `trim` to be `True`.

```
china.merge_layer(method = 'sum', weights = {
    'bins': 0.05, 'forest': 0.25, 'slope': 0.4, 'modis_forest': 0.3
}, trim = True)
```

### Loading shapes, extracting shape from mask

(WE NEED A DESCRIPTION OF WHAT THIS IS DOING)

Let us get province shapes from `cartopy` and save the path as `prov_path`. This can also be the path to user-supplied shape files.

```
prov_path = shpreader.natural_earth(resolution='10m', category='cultural', name = 'admin_1_states_provinces')
prov_path
```

Check attributes in the shapes contained in path `prov_path`

```
geodata.mask.shape_attribute(prov_path)
```

(DO WE NEED THIS PRINTOUT?)
Check out the get_shape() method:

```
geodata.mask.get_shape(
    path,
    key,
    targets=None,
    save_record=False,
    condition_key=None,
    condition_value=None,
    return_dict=False,
)
    Docstring:
    Take the path of the shape file, a attribute key name as the key to store in the output,
    and a list of targets such as provinces, return the shape of targets with attributes.

    While the targets list the exact names/value of key_name attributes of the shapes, users may
    also ignore it and use condition_key and condition_value to find the desired shapes.
    for example, the call below will find all the shapes of provinces that belongs to China:
        mask.get_shape(path_to_province_shapes, key = 'name',
            condition_key = 'admin', condition_value = 'China')

    We can also ignore condition, just take three provinces of China by naming them out.
    for example, the call below just take three provinces of China by naming them out:
        mask.get_shape(prov_path, key = 'name_en',
                    targets = ['Jiangsu', 'Zhejiang', 'Shanghai'])

    path (str): string path to the shapefile
    key_name (str): key name for the shapefile, can be checked through shape_attribute(path)
    targets (list): target names, if not specified, find all shapes contained in the shapefile.
    save_record (bool): if the records for the shape should be returned as well.
    condition_key (str): optional: the attribute as another condition
    condition_value (str)：optional: the value that the condition_key must equals to in the
        shape record.
    return_dict (str): if a dictionary will be returned, otherwise return a dataframe.
        False by default.

    return (dict or pandas.dataframe) the array of feature/shapes extracted from the shapefile
```
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
china_shapes
```

Extract the shapes from the merged_mask:

```
>>> china.extract_shapes(china_shapes, crop = True)
>>> china
Mask China:
4 layers: ['bins', 'forest', 'slope', 'modis_forest'] .
Merged_mask merged/flattened.
3 shape_mask: ['Zhejiang', 'Shanghai', 'Jiangsu'].
Mask has not been saved/updated.
```

## SHAPE AS LAYER (PERHAPS A MORE INTUITIVE NAME)

(USE ENVIRONMENTAL PROTECTED AREAS HERE INSTEAD OF PROVINCE)

This is different from shape extractions, as we will simply treat one shp file as a layer, instead of grabbing the merged mask within that shape.

The `add_shape_layer` method take in a dictionary of shapes, a resolution of the result raster with that shape.

```
china.add_shape_layer({'Jiangsu': china_shapes['Jiangsu']},
                      resolution = (5000, 5000))
show(china.layers['Jiangsu'])
```

The previous method call is great when we do not have any layers in the mask object at all. However, since we already have bins, forest, and slope... as the layers, we may want to make the shape layers the similar resolution, and similar boundary with the other layers. We can just specify a `reference_layer` in the method call. We can set it to a layer name in the mask, and the result layer resolution and bounds of the shape layer will be same with that specified layer in the mask object.

```
china.add_shape_layer(china_shapes,
                      reference_layer = 'bins')
show(china.layers['Jiangsu'])
```

## Saving the mask

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


## Loading a previously saved mask

```
>>> china_2 = geodata.mask.load_mask("china")
>>> china_2
Mask china:
7 layers: ['bins', 'forest', 'Jiangsu', 'modis_forest', 'Shanghai', 'slope', 'Zhejiang'] .
Merged_mask merged/flattened.
3 shape_mask: ['Jiangsu', 'Shanghai', 'Zhejiang'].
Mask has been saved.
```

### Possible errors to avoid (CAN WE PROVIDE THIS INFORMATION IN A DIFFERENT WAY)

If you create another object `china_2` that opens the raster object `china` is accessing, and then try to save the original `china` without using `china_2.close_files()`, you should expect an error because Python does not want you to rewrite a file that is used by another program. Therefore, `china_2.close_files()` or `china_2.save_mask(close_files = True)` make sures that only one mask is having access to the files. `close_files()` will close all the layers in china_2 and make that mask object un-savable. Therefore, it is best to avoid having multiple mask objects accessing the same files.

Please see the jupyter notebook for more details.
