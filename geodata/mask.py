## Copyright 2021 Jiahe Feng (Davidson Lab).

## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License as
## published by the Free Software Foundation; either version 3 of the
## License, or (at your option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.

## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.


"""
GEODATA
Geospatial Data Collection and "Pre-Analysis" Tools
"""

import numpy as np
import matplotlib.pyplot as plt
import cartopy.io.shapereader as shpreader
import logging
logger = logging.getLogger(__name__)
import geopandas as gpd
import shapely

import os, sys
import rasterio as ras
from rasterio.plot import plotting_extent
from rasterio.io import MemoryFile
from rasterio import mask as rmask
from rasterio.merge import merge

import xarray as xr
from . import config, datasets


class Mask(object):
    """
    A Module to create, manipulate, and load geodata mask object that takes geo tif or shp files as input.

    name (str): name of the mask object
    layers (dict): dictionary that stores the layers of the mask with 
            their names and ras.DatasetReader values
    merged_mask (ras.DatasetReader): processed mask after merging and flattening layers
    shape_mask (dict): dictionary that stores the extracted shape after
            masking the shape on the merged_mask (by default, or on certain layers if sepcified)
            with the shape names and ras.DatasetReader values
    saved (boolean): whether the mask has been saved/updated
    """

    ### INITIALIZATION

    def __init__(self, name, layer_path = None, layer_name = None, mask_dir = config.mask_dir, **kwargs):
        """
        Initialize a new mask object. Layer_path will take the file name without extension as the default layer name, but when
        layer name(s) is/are specified, it will take the layer_name(s) as the name for new added layer(s). Layer path can also
        be a dictionary of {layer_name(str): layer_path(str)} mapping.
        
        name (str): Name for the new mask object
        layer_path (str, list, or dictionary): Path(s) of new layers to be added to the mask object
        layer_name: (str, list): Names for Mask.layers created from the layer path
        mask_dir: (str): the path to where the mask object will be saved/stored. 
            By default, it should be the mask path in config.py

        """
        self.name = name
        self.layers = {}

        self.merged_mask = None 
        self.shape_mask = {}
        
        self.saved = False
        self.mask_dir = mask_dir
        
        if layer_path: #if there is layer path input
            self.add_layer(layer_path, layer_name, **kwargs)
            
    def __repr__(self):
        """Display representation of the mask object"""

        if not self.merged_mask: 
            mmask = 'No merged_mask ready. \n'
        else:
            mmask = 'Merged_mask merged/flattened. \n'

        if not self.shape_mask:
            smask = 'No shape has been extracted. \n'
        else:
            smask = f"{len(self.shape_mask)} shape_mask: {list(self.shape_mask.keys())}. \n"

        if not self.saved:
            is_save = "Mask has not been saved/updated. \n"
        else:
            is_save = "Mask has been saved. \n"
        
        return (
        f"Mask {self.name}: \n" +
        f"{len(self.layers)} layers: {list(self.layers.keys())}.\n" + mmask + smask 
        + is_save
        )

    def __str__(self):
        """Print/str representation of the object"""
        return (f"Mask {self.name}")

    ### LOADING LAYERS    

    @classmethod            
    def load_mask(self, name, mask_dir = config.mask_dir):
        """
        Load a previously saved mask object

        Parameters:
        name: name for the mask
        mask_dir: directory to look for previously saved mask file, and the path to where the mask object will be updated. 
            By default, it should be the mask path in config.py
        """
        obj_dir = os.path.join(mask_dir, name)

        #create new mask object
        prev_mask = Mask(name = os.path.basename(obj_dir))
        prev_mask.mask_dir = mask_dir
        
        #load the layers
        layer_path = os.path.join(obj_dir, 'layers')
        if (os.path.isdir(layer_path)):
            for i in os.listdir(layer_path):
                prev_mask.layers[i.split('.')[0]] = ras.open(os.path.join(layer_path, i))
                logger.info(f"Layer {i.split('.')[0]} loaded to the mask {prev_mask.name}.")
        else:
            logger.info("No previously saved mask found for mask %s", name)

        #load the merged layer, if there is one     
        merge_path = os.path.join(obj_dir, 'merged_mask')
        merge_file_path = os.path.join(merge_path, 'merged_mask.tif')
        if (os.path.exists(merge_file_path)):
            prev_mask.merged_mask = ras.open(merge_file_path)
            logger.info(f"Merged_mask loaded to the mask {prev_mask.name}.")
        else:
            logger.info("No Merged Mask found.")

        #load the shape masks  
        shape_path = os.path.join(obj_dir, 'shape_mask')
        if (os.path.isdir(shape_path)):
            for i in os.listdir(shape_path):
                prev_mask.shape_mask[i.split('.')[0]] = ras.open(os.path.join(shape_path, i))
                logger.info(f"Shape mask {i.split('.')[0]} loaded to the mask {prev_mask.name}.")
        
        prev_mask.saved = True

        return prev_mask
    
    @classmethod
    def open_tif(self, layer_path, show = True, close = False):
        """
        Openning a layer using rasterio given a file path and store it as openning_layer attribute. 
        This method does not add the layer to the mask object.

        layer_path (str): the path to the tif file
        show (bool): if the method will plot the raster
        close (bool): if the method will close the raster afterward, prevent possible permission errors.
        """
        openning_tif = ras.open(layer_path)
        if show: Mask.show(openning_tif)
        if not close: 
            logger.info("Please remember to close the file with .close()")
        else:
            openning_tif.close()
        return openning_tif

    def _add_layer(self, layer_path, layer_name = None, replace = True, trim_raster = False):
        """
        Add a layer to the mask, this method incorporate CRS convertion. 

        layer_path (str): the path to the tif file
        layer_name (str): the layer name, by default None and the layer name would be the file name
        replace (bool): if the layer with same name will replace the old one
        trim_raster (bool): if the method should trim the all-empty row/column border of the raster

        """
   
        if not layer_name:
            layer_name = os.path.basename(layer_path).split(".")[0]
            
        #replace layer by default
        if layer_name in self.layers:
            if replace == True:
                self.layers[layer_name].close()
                del self.layers[layer_name] #delete old layer from memory
                logger.info(f"Overwriting existing layer {layer_name}.")
            else:
                raise ValueError("Layer name %s already existed in mask %s, "
                           "please change the name or replace the existing one with replace = True", 
                           layer_name, self.name)
        
        new_raster = ras.open(layer_path)
        
        #check if CRS is lat-lon system
        if ras.crs.CRS.from_string('EPSG:4326') != new_raster.crs:
            self.layers[layer_name] = Mask.reproject_raster(new_raster, trim_raster = trim_raster)
        else: 
            self.layers[layer_name] = new_raster
        
        logger.info(f"Layer {layer_name} added to the mask {self.name}.")
        self.saved = False
        
    def add_layer(self, layer_path, layer_name = None, **kwargs):
        """Add a layer to the mask from given path(s). Called the _add_layer() helper method"""

        if isinstance(layer_path, dict): #if layer_path input is dictionary
            for name, val in layer_path.items():
                self._add_layer(val, layer_name = name, **kwargs)
            return

        if isinstance(layer_path, str): #if layer path is string
            layer_path = [layer_path]
            layer_name = [layer_name]

        for i in range(len(layer_path)): #add the layers, using the _add_layer() method
            self._add_layer(layer_path[i], layer_name[i], **kwargs)

        return

    def remove_layer(self, name):
        """Remove a layer in the mask given name."""
        if name in self.layers:
            self.layers[name].close()
            del self.layers[name]
        else:
            raise KeyError("No layer name %s found in the mask.", name)
        self.saved = False
        

    def rename_layer(self, layer_name, new_name):
        """Rename a layer in the mask given name."""
        if layer_name in self.layers:
            self.layers[new_name] = self.layers[layer_name]
            del self.layers[layer_name]
            self.saved = False
        else:
            raise KeyError("No layer name %s found in the mask.", layer_name)


    ### USEFUL HELPER METHODS FOR LAYERS

    def get_res(self, layers = None, product = False):
        """
        Get resolutions of layers.

        layers (list): the list of layers to get resolution. If not specified, select all layers in the object.
        product (bool): if the user want to see the product of the grid cell height and width.
        """
        res = {}
        
        if not layers: 
            layers = self.layers
        else:
            layers = {k: self.layers[k] for k in layers}
            
        for k, val in layers.items():
            res[k] = val.res

        if product: #return grid cell area (not geographical)
            for k in layers.keys():
                res[k] = res[k][0] * res[k][1]

        return res
    
    def get_bounds(self, layers = None):
        """Get boundary information of layers"""
        bounds = {}
    
        if not layers: 
            layers = self.layers
        else:
            layers = {k: self.layers[k] for k in layers}
        
        for k, val in layers.items():
            bounds[k] = val.bounds
        return bounds
    
    def find_inner_bound(self, layers = None):
        """Find the bound where all layers overlap"""

        layer_bound = self.get_bounds(layers = layers)
        min_left, min_bottom = (max([layer_bound[b][num] for b in layer_bound]) for num in [0, 1])
        min_right, min_top = (min([layer_bound[b][num] for b in layer_bound]) for num in [2, 3])

        return (min_left, min_bottom, min_right, min_top)
    
    @classmethod
    def create_temp_tif(self, arr, transform, open_raster = True, **kwargs):
        """
        Create a ras.DatasetReader object openning a temporary rasterio file
        Parameters:
        arr (np.array): np.array contains values for that layer
        transform (affine.Affine): affine transform information for that layer
            arr, transform input can be the output for ras.merge.merge and ras.mask.mask)
        open_raster: if the raster will be opened
        """
        with MemoryFile() as memfile: #using memory file to avoid creating new files
            with ras.open(
                memfile.name,
                'w',
                driver='GTiff',
                height=arr.shape[0],
                width=arr.shape[1],
                count=1,
                dtype=arr.dtype,
                compress='lzw',
                crs='+proj=latlong',
                transform=transform,
            ) as dst:
                dst.write(arr, 1)
                
            if open_raster == True:
                return ras.open(memfile.name)
            else:
                return memfile.name
            
    @classmethod     
    def _save_opened_raster(self, raster, path):
        """Helper method to close the opened raster and save it to the given path"""
        arr, transform = raster.read(1), raster.transform
        raster.close()
        Mask.save_raster(arr, transform, path)


    @classmethod
    def save_raster(self, arr, transform, path, **kwargs):
        """
        Given np.array and Affine transform, save the raster as tif file to the path
        https://rasterio.readthedocs.io/en/latest/topics/writing.html

        arr (np.array): np.array contains values for that layer.
        transform (affine.Affine): affine transform information for that layer
            arr, transform input can be the output for ras.merge.merge and ras.mask.mask)
        path (str): the path to where the tif file is saved.
        """
        with ras.open(
            path,
            'w',
            driver='GTiff',
            height=arr.shape[0],
            width=arr.shape[1],
            count=1,
            dtype=arr.dtype,
            compress='lzw',
            crs='+proj=latlong',
            transform=transform,
        ) as dst:
            dst.write(arr, 1)

    @classmethod
    def crop_raster(self, raster, bounds, lat_lon_bounds = True):
        """
        Crop a raster given geographical coordinates or its array index bounds

        raster (ras.DatasetReader): raster openner / opened rasterio tif dataset
        bounds (tuple: (left, bottom, right, top)): geographical coordinates or array index bounds
        lat_lon_bounds (boolean): if True (by default), the bounds tuple are latitude-longitude coords input
        """

        if lat_lon_bounds: #if bounds are lat/long values
            window = ras.windows.from_bounds(bounds[0], bounds[1], bounds[2], bounds[3], raster.transform)
        else: #bounds are index
            window = ras.windows.Window.from_slices((bounds[0], bounds[1]), (bounds[2], bounds[3]))

        with MemoryFile() as memfile: #using memory file to avoid creating new files

            kwargs = raster.meta.copy()
            kwargs.update({
                'height': window.height,
                'width': window.width,
                'transform': ras.windows.transform(window, raster.transform)})

            with ras.open(memfile.name, 'w', compress='lzw', **kwargs) as dst:
                dst.write(raster.read(window=window))
                dst.close()

            return ras.open(memfile.name)
    
    def crop_layer(self, layer_name, bounds, lat_lon_bounds = True, dest_layer_name = None):
        """Crop a layer of mask object using Mask.crop_raster() method"""
        if not dest_layer_name: dest_layer_name = layer_name
        self.layers[dest_layer_name] = Mask.crop_raster(self.layers[layer_name], bounds, lat_lon_bounds)
        self.saved = False

    @classmethod
    def trim_raster(self, raster):
        """
        Remove the all-zero columns and rows at the border of the raster and returns the trimmed raster.
        Do not remove the all-zero cols or rows in the middle of the valid values.

        for example: if the raster.read(1) (the array values) is the following:
            (np.array([[0, 0, 9, 0, 0, 0, 0],
                    [0, 0, 1, 2, 0, 0, 0],
                    [0, 0, 2, 3, 4, 0, 0],
                    [0, 0, 0, 0, 0, 0, 0],
                    [0, 0, 1, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0, 0]]))
        The bounds of valid values will be (2, 4, 0, 4);
        The all zero columns and rows at the border of the array will be removed.
        """

        arr = raster.read(1)
        
        all_zero_col = np.argwhere(np.all(arr[..., :] == 0, axis=0)).reshape(-1)
        all_zero_row = np.argwhere(np.all(arr[:,...] == 0, axis=1)).reshape(-1)

        l_idx, t_idx = 0, 0
        r_idx, b_idx = arr.shape[1] - 1, arr.shape[0] - 1
        
        #find left and right index where all empty column start and ends
        if len(all_zero_col) != 0:
            if all_zero_col[0] == 0:
                for i in range(len(all_zero_col) - 1):
                    l_idx += 1
                    if all_zero_col[i + 1] != all_zero_col[i] + 1:
                        break
            if all_zero_col[-1] == r_idx:
                for i in range(len(all_zero_col) - 1, 0, -1):
                    r_idx -= 1
                    if all_zero_col[i] != all_zero_col[i - 1] + 1:
                        break
        
        #find top and bottom index where all empty column start and ends
        if len(all_zero_row) != 0:
            if all_zero_row[0] == 0:
                for i in range(len(all_zero_row) - 1):
                    t_idx += 1
                    if all_zero_row[i + 1] != all_zero_row[i] + 1:
                        break

            if all_zero_row[-1] == b_idx:
                for i in range(len(all_zero_row) - 1, 0, -1):
                    b_idx -= 1
                    if all_zero_row[i] != all_zero_row[i - 1] + 1:
                        break

        bounds = (t_idx, b_idx, l_idx, r_idx)

        return Mask.crop_raster(raster, bounds, lat_lon_bounds = False)

    def trim_layer(self, layer_name, dest_layer_name = None):
        """Trim a layer of mask object using Mask.trim_raster() method"""
        if not dest_layer_name: dest_layer_name = layer_name
        self.layers[dest_layer_name] = Mask.trim_raster(self.layers[layer_name])
        self.saved = False

    @classmethod
    def reproject_raster(self, src, dst_crs = 'EPSG:4326', trim_raster = False, **kwargs):
        """
        convert non-lat-lon crs tif file to lat-lon crs

        src (ras.DatasetReader): the source raster
        dst_crs (str): by default, we want to make the destination raster lat-lon
        trim_raster (bool): False by default. If True, we will trim the empty border to save space. 
        """

        transform, width, height = ras.warp.calculate_default_transform(
            src.crs, dst_crs, src.width, src.height, *src.bounds)
        kwargs = src.meta.copy()
        kwargs.update({
            'crs': dst_crs,
            'transform': transform,
            'width': width,
            'height': height
        })

        #write it to another file: the CRS corrected one
        #rasterio.readthedocs.io/en/latest/topics/reproject.html
        with MemoryFile() as memfile:
            with ras.open(memfile.name, 'w', compress='lzw', **kwargs) as dst:
                for i in range(1, src.count + 1):
                    ras.warp.reproject(
                        source = ras.band(src, i),
                        destination = ras.band(dst, i),
                        src_transform = src.transform,
                        src_crs = src.crs,
                        dst_transform = transform,
                        dst_crs = dst_crs,
                        resampling = ras.warp.Resampling.nearest)

            logger.info(f"Raster {src.name} has been reprojected to lat-lon CRS.")    
            return_ras =  ras.open(memfile.name)
            
            if trim_raster:
                return Mask.trim_raster(return_ras)  

            return return_ras

    @classmethod          
    def binarize_raster(self, raster, values):
        """
        Map values to one in a new binary raster layer

        raster (ras.DatasetReader): the source raster
        values (list): the list of numberic values in the raster array to filter out as 1
        
        """
        new_arr = np.isin(raster.read(1), values) * 1
        return Mask.create_temp_tif(new_arr.astype(np.int8), raster.transform)
     
    def binarize_layer(self, layer_name, values, dest_layer_name = None):
        """binarize a layer of mask object using Mask.binarize_raster() method"""
        if not dest_layer_name: dest_layer_name = layer_name
        self.layers[dest_layer_name] = Mask.binarize_raster(self.layers[layer_name], values)
        self.saved = False

    ### MERGING AND FLATTENING
    
    def _sum_method(merged_data, new_data, merged_mask, new_mask, **kwargs):
        """The sum method will add up the values from all the layers. We can also customize the weights. 
        Behind the scene of the sum method: it multiplys each layers with the corresponding weight, 
        and add the in-memory temporary layers together"""
        
        if len(np.unique(merged_data)) == 1:
            mask = np.ones(merged_data.shape, dtype=bool)
            np.add(np.zeros(merged_data.shape), new_data.data, out=merged_data, where = mask, casting="unsafe")
        else:
            np.add(merged_data, new_data.data, out=merged_data, casting="unsafe")
        
    def _and_method(merged_data, new_data, merged_mask, new_mask, **kwargs):
        """By default, the merge_layer method will use a binary and method: 
        if any of the 4 grid cells of the 4 layers at the same location have 0, 
        then the returned self.merged_layer will also have 0 at that location. 
        Therefore, if all the layers indicate that a land is not unavailable (!=0), 
        the merged result will have value 1. """

        if len(np.unique(merged_data)) == 1:
            mask = np.ones(merged_data.shape, dtype=bool)
        else:
            mask = merged_data != 0
        
        np.copyto(merged_data, new_data.data, where=mask, casting="unsafe")
        
    def merge_layer(self, method = 'and', weights = None, 
                    layers = None, trim = False, show = True, reference_layer = None, 
                    attribute_save = True, **kwargs):
        """
        Merge multiple layers together and flatten it using either 'and' and 'sum' method. By default, 'and' method
        is used, and we would save the result to self.merged_mask attribute.

        Adjusted from rasterio's documentation on merging: https://rasterio.readthedocs.io/en/latest/api/rasterio.merge.html
        Geospatial bounds and resolution of a new output file in the units of the input file coordinate reference system
        may be provided, but by default, we will use the layer with the best resolution for the output bounds/resolution, 
        unless a reference layer is provided.

        method (str): either "sum" or "and" method, the difference is explained in the docstring for _sum_method and _and_method
        weights (dict): when the "and" method is used, the weights for each layer. None by default, and in this case the weight will
            be 1 for all layers.
        layers (list): a list of layers to be merged. None by default and the method will automatically perform merging on all the layers
            contained in the same mask object.
        trim (bool): if we will trim the area where all the layers must overlap, False by default.
        show (bool): whether the merged_mask output will be plotted, True by default.
        reference_layer (str): the name of the layer which we use as the output bounds/resolution.
        attribute_save (bool): whether the method save the output to self.merged_mask, True by default.
        **kwargs: other parameters for rasterio.merge.merge(). Use the link above to find all possible parameters.

        """
        
        if not layers: layers = list(self.layers.keys())
        
        #specify layers to be used in the merged process
        layers = {k: self.layers[k] for k in layers}
        
        if method == 'and':
            arr, aff = merge(list(layers.values()), method = Mask._and_method, **kwargs)
            
        elif method == 'sum':
            if not weights:
                logger.warning('No weight dictionary provided, all the layers for merging will have weights of 1 by default')
                weights = {k: 1 for k in layers.keys()}
            else:
                for k in layers.keys():
                    if k not in set(weights.keys()):
                        weights[k] = 1
            
            if not reference_layer: #if no reference layer is provided, using the layer with best resolution

                #find all resolutions
                resolutions = {}
                for k, val in layers.items():
                    resolutions[k] = val.res

                #find the layer with best resolution   
                best_res = float('inf')
                for k, res in resolutions.items():
                    if res[0] * res[1] < best_res:
                        reference_layer = k
                        best_res = res[0] * res[1]

                temp_layers = {}
                for lay_name, lay_weight in weights.items():
                    temp_layers[lay_name] = Mask.create_temp_tif(layers[lay_name].read(1) * lay_weight,
                                                    layers[lay_name].transform)   

            #make sure highest resolution layer is the first input for ras.merge
            merged_lst = [temp_layers[reference_layer]]
            temp_layers.pop(reference_layer)
            merged_lst += list(temp_layers.values())
        
            arr, aff = merge(merged_lst, method = Mask._sum_method, **kwargs)
            [r.close() for r in merged_lst]
            
        return_ras = Mask.create_temp_tif(arr[0], aff)
        
        if trim:
            bounds = self.find_inner_bound()
            return_ras = Mask.crop_raster(return_ras, bounds)

        if show: 
            Mask.show(return_ras, title = 'Merged Mask')

        if attribute_save == True:
            if self.merged_mask:
                logger.info("Overwriting current merged_mask.")
            self.merged_mask = return_ras
            logger.info("Merged Mask saved as attribute 'merged_mask'.")
            
        else:
            return return_ras

    def clear_merge_layer(self):
        self.merged_mask == None
        
    ### LOADING SHAPES, CREATING SHAPE MASK, EXTRACTING SHAPE FROM MASK
    
    def shape_attribute(path):
        """
        Get the first shape objects (key value pairs) in that path to check attributes
        Knowing the key is important as we need to specify key names in get_shapes()

        path (str): string path to the shapefile
        """
        return next(shpreader.Reader(path).records()).attributes 
    
    def get_shape(path, key, targets = None, save_record = False, 
                  condition_key = None, condition_value = None,
                  return_dict = False):
        """
        Take the path of the shape file, a attribute key name as the key to store in the output, 
        and a list of targets such as provinces, return the shape of targets with attributes.

        While the targets list the exact names/value of key_name attributes of the shapes, users may
        also ignore it and use condition_key and condition_value to find the desired shapes.
        for example, the call below will find all the shapes of provinces that belongs to China:
            Mask.get_shape(path_to_province_shapes, key = 'name', 
                condition_key = 'admin', condition_value = 'China')
            
        We can also ignore condition, just take three provinces of China by naming them out.
        for example, the call below just take three provinces of China by naming them out:
            Mask.get_shape(prov_path, key = 'name_en', 
                         targets = ['Jiangsu', 'Zhejiang', 'Shanghai'])

        path (str): string path to the shapefile
        key_name (str): key name for the shapefile, can be checked through shape_attribute(path)
        targets (list): target names, if not specified, find all shapes contained in the shapefile.
        save_record (bool): if the records for the shape should be returned as well.
        condition_key (str): optional: the attribute as another condition
        condition_value (str)：optional: the value that the condition_key must equals to in the shape record.
        return_dict (str): if a dictionary will be returned, otherwise return a dataframe. False by default.
        """
        reader = shpreader.Reader(path)

        shapes = []
        shape_keys = []
        records = {}

        if (condition_key and not condition_value) or (condition_key and not condition_value):
            raise ValueError("Condition key and condition value must be both null, or have values.")

        #if not conditional key and value are specified
        if (not condition_key and not condition_value):
            if not targets:
                for r in reader.records():
                    r_key = r.attributes[key]
                    shape_keys.append(r_key)
                    shapes.append(r.geometry)
                    records[r_key] = r.attributes
            else: #target list provided
                for r in reader.records():
                    r_key = r.attributes[key]
                    if r_key in targets:
                        shape_keys.append(r_key)
                        shapes.append(r.geometry)
                        records[r_key] = r.attributes

        else:
            #with extra conditions
            if not targets:
                for r in reader.records():
                    r_key = r.attributes[key]
                    if r.attributes[condition_key] == condition_value:
                        shape_keys.append(r_key)
                        shapes.append(r.geometry)
                        records[r_key] = r.attributes
            else: #target list provided
                for r in reader.records():
                    r_key = r.attributes[key]
                    if r.attributes[condition_key] == condition_value and r_key in targets:
                        shape_keys.append(r_key)
                        shapes.append(r.geometry)
                        records[r_key] = r.attributes

        shapes = gpd.GeoDataFrame({key: shape_keys, 'shapes': shapes})
        
        if return_dict: shapes = shapes.set_index(key).to_dict()['shapes']

        if save_record:
            return shapes, records

        return shapes
    
    def add_shape_layer(self, shapes, reference_layer = None, resolution = None, 
                        combine_shape = False, combine_name = None, invert = True, **kwargs):
        """
        Add shapes to the mask layers. This is different from shape extractions, 
        as we will simply treat one shp file as a layer, instead of grabbing the merged mask within that shape.
        This method take in a dictionary of shapes, a resolution of the result raster with that shape, and add
        the shape to the mask object. We can also use a reference layer that is present in the mask object to avoid
        manuelly finding resolution.

        shapes (dict): a dictionary of shapes
        reference_layer (str): name to the layer which we use its bounds/resolution
        resolution (tuple): tuple of (width, height). If specified with a reference layer, we ignore the resolution
            of the referenced layer and use the input resolution.
        combine_shape (bool): if the user want to combine the shapes as one shape. Only one layer will be added as a result.
        combine_name (str)： the name of the combined shape if combine_shape is True
        invert (bool): True if we want to have 1 inside of shape and 0 outside of the shape in the added layer. True by default.
        """

        if reference_layer:
            shape_transform = self.layers[reference_layer].transform
            if not resolution:
                resolution = self.layers[reference_layer].shape
            else:
                logger.info(f"The resolution of the new layers(s) will be {resolution}.")

        else:
            if not resolution:
                raise TypeError("lPlease provide the tuple of (width, height) for the resolution of the new shape raster.")

        if combine_shape:
            if not combine_name: raise ValueError("Please specify combined shape name.")
                
            shapes = {combine_name: shapely.ops.unary_union(shapes.values())}

        for name, shp in shapes.items():
            if not isinstance(shp, shapely.geometry.multipolygon.MultiPolygon):
                #make sure that the shape is a geometry.multipolygon
                shp = shapely.geometry.multipolygon.MultiPolygon([shp])
            shape_bounds = ras.features.bounds(shp)
            if not reference_layer:  
                shape_transform = ras.transform.from_bounds(*shape_bounds, *resolution)
            arr = ras.features.geometry_mask(shp, resolution, shape_transform, invert = invert, **kwargs) * 1
            self.layers[name] = Mask.create_temp_tif(arr.astype(np.int8), shape_transform)
            logger.info(f"Layer {name} added to the mask {self.name}.")


    def extract_shapes(self, shapes, layer = None, 
                       combine_shape = False, combine_name = None,
                       show = False, attribute_save = True, **kwargs):
        """
        Extract the shapes on the (by default) merged_mask layer, and save the result shape raster in 'shape_mask' attribute
        
        shapes (dict): a dictionary of shapes
        layer (str): the layer to extract shape from. By default the layer should be None and we want to extract the shapes from
            the shape_mask layer.
        combine_shape (bool): if the user want to combine the shapes as one shape. Only one layer will be added as a result.
        combine_name (str)： the name of the combined shape if combine_shape is True
        show (bool): if we want to plot the result shapes. True by default.
        attribute_save (bool): if we want to save the shapes to the shape_mask attributes. True by default.
        """

        if not layer:
            if not self.merged_mask:
                raise ValueError("No merged_mask ready for the mask object. Please create one using merge_layer() call"
                              ", or use a layer input by sepcifiying layer = [layer's name]")
            else:
                layer = self.merged_mask

        else:
            layer = self.layers[layer]

        if isinstance(shapes, dict) == False: 
            raise ValueError("shapes input must be dictionary (key-value pair).")

        if combine_shape:
            if not combine_name: raise ValueError("Please specify combined shape name.")
                
            shapes = {combine_name: shapely.ops.unary_union(shapes.values())}

        return_shape = {}

        #loop through shapes to extract and create temp tif raster
        for name, shp in shapes.items():
            if isinstance(shp, list) == False: shp = [shp]
            arr, aff = rmask.mask(layer, shp, nodata = 0, **kwargs)
            raster = Mask.create_temp_tif(arr[0], aff)
            return_shape[name] = raster
            if attribute_save:
                if name in self.shape_mask:
                    self.shape_mask[name].close()
                    logger.info(f"[Overwritten] Extracted shape {name} added to attribute 'shape_mask'.") 
                else:
                    logger.info(f"Extracted shape {name} added to attribute 'shape_mask'.") 

        if show:
            [Mask.show(r, title=name) for name, r in return_shape.items()]

        if attribute_save:
            self.shape_mask.update(return_shape)
            self.saved = False
 
        else: 
            return return_shape

    def remove_shapes(self, names):
        """Remove a shape from shape_mask dictionary"""
        for n in names:
            if n not in self.shape_mask.values():
                raise KeyError("Shape mask %s not found in the object.", n)
            self.shape_mask[n].close()
            del self.shape_mask[n]
                
    def load_shape_xr(self, names = None):
        """Using xarray to load the shape masks."""
        if self.saved == False:
            raise ValueError(f"The Mask object has not been saved, please call save_mask() first.")
        shape_path = os.path.join(self.mask_dir, self.name, 'shape_mask')
        xrs = {}
        if not names: names = self.shape_mask.keys()
        for n in names:
            xrs[n] = xr.open_rasterio(os.path.join(shape_path, (n + '.tif')))
        return xrs

    ### SAVING MASK

    def close_files(self):
        """Close all the opened rasters. This method will disable further save_mask() call."""
        [i.close() for i in self.layers.values()]
        if self.merged_mask: self.merged_mask.close()
        if self.shape_mask:
            [i.close() for i in self.shape_mask.values()]

    def save_mask(self, name = None, mask_dir = None, close_files = False):
        """
        Save a mask object to the directory. By default, the directory should be the mask_dir in config.py.
        It is recommended that multiple masks would not be opened at the same time to avoid permission issues.
        
        name (str): name of the mask object to be saved.
        mask_dir (path): path to save the mask to. By default, use the value in mask_dir attribute.
        close_files (bool): if close all the opened raster after saving each. This will disable further save_mask() call.
            by default it is False.
        """
    
        #if user want to use another name to save the mask
        if not name: name = self.name
        extension = '.tif'
        
        if not mask_dir: mask_dir = self.mask_dir

        # Ensure mask directory exists
        if not os.path.isdir(mask_dir):
            os.mkdir(mask_dir)

        #get mask object path
        obj_path = os.path.join(mask_dir, self.name)
        if not (os.path.isdir(obj_path)):
            os.mkdir(obj_path)

        #layer path
        layer_path = os.path.join(obj_path, 'layers')
        if not (os.path.isdir(layer_path)):
            os.mkdir(layer_path)
        else: #if directory do not exist, remove all the previously saved content
            for f in os.listdir(layer_path):
                if f.split('.')[0] not in self.layers.keys():
                    os.remove(os.path.join(layer_path, f))
        #update layer
        if self.layers:
            for name, raster in self.layers.items():
                Mask._save_opened_raster(raster, os.path.join(layer_path, name + extension))
                if not close_files: self.layers[name] = ras.open(os.path.join(layer_path, name + extension))

        #merge_layer path
        merge_path = os.path.join(obj_path, 'merged_mask')
        merge_tif_path = os.path.join(merge_path, 'merged_mask.tif')

        if self.merged_mask:
            if not(os.path.isdir(merge_path)):
                os.mkdir(merge_path)
            Mask._save_opened_raster(self.merged_mask, merge_tif_path)
            if not close_files: self.merged_mask = ras.open(merge_tif_path)
        else: 
            if os.path.isfile(merge_tif_path):#deleted previously one
                os.remove(merge_tif_path)

        #shape path
        shape_path = os.path.join(obj_path, 'shape_mask')
        if not(os.path.isdir(shape_path)):
            os.mkdir(shape_path)
        else:
            for f in os.listdir(shape_path):
                if f.split('.')[0] not in self.shape_mask.keys():
                    os.remove(os.path.join(shape_path, f))
        if self.shape_mask:
            for name, raster in self.shape_mask.items():
                Mask._save_opened_raster(raster, os.path.join(shape_path, name + extension))
                if not close_files: self.shape_mask[name] = ras.open(os.path.join(shape_path, name + extension))
            
        self.saved = True
        logger.info(f"Mask {self.name} successfully saved at {self.mask_dir}")
        
    ### VISUALIZATION METHOD
    @classmethod
    def show(self, raster, title=None, colorbar = True, **kwargs):
        """
        Plot a rasterio file given ras.DatasetReader

        raster (ras.DatasetReader): rasterio file opener, the value in the layers, shape_mask, and merged_mask attribute.
        title (str): the title of the plot
        colorbar (bool): if show the legend for the values of the plot. True by default
        """
        plt.imshow(raster.read(1), 
                extent = plotting_extent(raster.read(1), raster.transform),
                **kwargs)
        if title == True: title = raster.name
        plt.title(title)
        if colorbar: plt.colorbar()
        plt.show()
    
    @classmethod
    def show_all(self, r_dict, **kwargs):
        """
        show all the layers given a dictionary
        Example use case: Mask.show_all(my_mask.layers); Mask.show_all(my_mask.shape_mask)
        
        """
        [Mask.show(r, title=name, **kwargs) for name, r in r_dict.items()]
    
