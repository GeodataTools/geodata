# Mask Troubleshooting

Possible errors for the mask module and troubleshooting information.


### Object Attribute Property

Not important for its normal usage, but it may be helpful:

Although we have openning layers from direct input files in Rasterio, which only reads/opens/writes data on disk, sometimes we may have in-memory tif files in the `layers` attributes, and the product `merged_mask` and `shape_mask` are also in-memory temporary files. This is because we might create a new layer after automatic CRS conversion if we detect that the input file is not in latitude-longitude CRS, and in rasterio, a CRS conversion, merging-flattening, shapes on raster, cropping, and many other methods that make changes to the raster will requires creating a new file on disk, but we can avoid creating too many temporarily files and deleting them later by using memory files. Read more here: https://rasterio.readthedocs.io/en/latest/topics/memory-files.html


### Saving the mask

Recall that if you want to save a mask called `china`, you will call the following method.

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

### Loading a previously saved mask

```
>>> china_2 = geodata.mask.load_mask("china")
>>> china_2
Mask china:
7 layers: ['bins', 'forest', 'Jiangsu', 'modis_forest', 'Shanghai', 'slope', 'Zhejiang'] .
Merged_mask merged/flattened.
3 shape_mask: ['Jiangsu', 'Shanghai', 'Zhejiang'].
Mask has been saved.
```

### Possible errors to avoid

If you create another object `china_2` that opens the raster object `china` is accessing, and then try to save the original `china` without using `china_2.close_files()`, you should expect an error because Python does not want you to rewrite a file that is used by another program. Therefore, `china_2.close_files()` or `china_2.save_mask(close_files = True)` make sures that only one mask is having access to the files. `close_files()` will close all the layers in china_2 and make that mask object un-savable. Therefore, it is best to avoid having multiple mask objects accessing the same files.

Please see the jupyter notebook for more details.
