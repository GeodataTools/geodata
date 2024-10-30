# Copyright 2021 Jiahe Feng (Davidson Lab).
# Copyright 2023 Xiqiang Liu.

# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License as
# published by the Free Software Foundation; either version 3 of the
# License, or (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

"""Mask module of the geodata package. This module contains the Mask class, which
is used to create, manipulate, and load geodata mask object that takes GeoTIFF and shapefiles
as input. Additionally, this module also contains a set of functions that are used to
interact with raw raster data and shapefiles.
"""

import logging
import os
from typing import Iterable, Literal, Optional

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pyproj
import rasterio as ras
import shapely
import xarray as xr
from rasterio import mask as rmask
from rasterio.io import MemoryFile
from rasterio.merge import merge
from rasterio.plot import plotting_extent
from rioxarray import open_rasterio

from .config import MASK_DIR

logger = logging.getLogger(__name__)

try:
    from numba import njit, prange
except ImportError:
    from .utils import dummy_njit as njit

    prange = range
    logging.info("Numba not installed, custom functions on arrays are not vectorized.")

# NOTE: Shapely's had an API change that results in MultiPolygon
# being non-iterable. Use shp.geoms in those cases.
SHAPELY_NEW_API = tuple(int(i) for i in shapely.__version__.split(".")) > (1, 8, 0)


class Mask:
    """A class to create, manipulate, and load geodata mask object that takes GeoTIFF
    or shp files as input.

    Creating a new mask object. Layer_path will take the file name(s) without extension
    as the default layer name. If layer name(s) is/are specified, it will take the corresponding
    layer_name(s) as the name for new layer(s). Layer path can also be a dictionary of
    `{layer_name: layer_path(str)}` key-value pairs.

    Args:
        name (str): Name for the new mask object.
        layer_path (str or dict): Path(s) of new layers to beadded to the mask object.
            It can be None. In that case, nothing will be added.
        layer_name (str or list): Names for the incoming layers.
            If layer_path is a dictionary, layer_name will be ignored.
        mask_dir (str): The path to where the mask object would be saved/stored.
            By default, it should be the MASK_PATH constant in config.py.
        **kwargs: Additional arguments passed to the `add_layer` method, only used if layers need to be added.
    """

    MergingMethods = Literal["and", "sum"]

    def __init__(
        self,
        name: str,
        layer_path: Optional[str | list | dict] = None,
        layer_name: str | list | None = None,
        mask_dir: str = MASK_DIR,
        **kwargs,
    ):
        self.name = name
        self.layers: dict[str, ras.DatasetReader] = {}

        self.merged_mask = None
        self.shape_mask = {}

        self.saved = False
        self.mask_dir = mask_dir

        if layer_path:  # if there is layer path input
            self.add_layer(layer_path, layer_name, **kwargs)

    def __repr__(self):
        """Display representation of the mask object"""

        if not self.merged_mask:
            mmask = "No merged_mask ready. \n"
        else:
            mmask = "Merged_mask merged/flattened. \n"

        if not self.shape_mask:
            smask = "No shape has been extracted. \n"
        else:
            smask = (
                f"{len(self.shape_mask)} shape_mask: {list(self.shape_mask.keys())}. \n"
            )

        if not self.saved:
            is_save = "Mask has not been saved/updated. \n"
        else:
            is_save = "Mask has been saved. \n"

        return (
            f"Mask {self.name}: \n"
            + f"{len(self.layers)} layers: {list(self.layers.keys())}.\n"
            + mmask
            + smask
            + is_save
        )

    def _add_layer(
        self,
        layer_path: str,
        layer_name: Optional[str] = None,
        replace: bool = True,
        trim: bool = True,
        src_crs: Optional[str] = None,
        dest_crs: str = "EPSG:4326",
    ):
        if not os.path.exists(layer_path):
            raise FileNotFoundError(f"{layer_path} not found.")

        if not layer_name:
            layer_name = os.path.basename(layer_path).split(".")[0]

        # replace layer by default
        if layer_name in self.layers:
            if replace is True:
                self.layers[layer_name].close()
                del self.layers[layer_name]  # delete old layer from memory
                logger.info("Overwriting existing layer %s.", layer_name)
            else:
                raise ValueError(
                    f"Layer name {layer_name} already existed in mask {self.name}, please change the name or "
                    "replace the existing one with replace = True."
                )

        new_raster = ras.open(layer_path, "r+")

        # make sure that nodata value is 0
        new_raster.nodata = 0

        if not src_crs:
            src_crs = new_raster.crs

        if src_crs != dest_crs:
            # check if CRS is lat-lon system
            if ras.crs.CRS.from_string(dest_crs) != new_raster.crs:
                new_raster = reproject_raster(
                    new_raster, src_crs=src_crs, dst_crs=dest_crs, trim=trim
                )

        self.layers[layer_name] = new_raster

        logger.info("Layer %s added to the mask %s.", layer_name, self.name)
        self.saved = False

    def add_layer(
        self,
        layer_path: str | list | dict[str, str],
        layer_name: str | list | None = None,
        replace: bool = True,
        trim: bool = True,
        src_crs: Optional[str] = None,
        dest_crs: str = "EPSG:4326",
    ):
        """Add a layer to the mask from given path(s) by calling the _add_layer() method.

        Args:
            layer_path (str or list or dict): The path(s) to the layer file(s)
            layer_name (str or list): The layer name(s), if names weren't specificd in
                `layer_path`.
            replace (bool): Whether the layer input with same name will replace the old one
                By default True.
            trim (bool): Whether the method will trim the all-empty row/column border
                of the raster. By default, this is true.
            src_crs (str): The source file CRS
            dest_crs (str): The destination CRS.
                By default it is 'EPSG:4326' lat lon coordinate system.
        """

        if isinstance(layer_path, dict):  # if layer_path input is dictionary
            for name, val in layer_path.items():
                self._add_layer(
                    val,
                    layer_name=name,
                    replace=replace,
                    trim=trim,
                    src_crs=src_crs,
                    dest_crs=dest_crs,
                )
            return

        if isinstance(layer_path, str):  # if layer path is string
            layer_path = [layer_path]
            layer_name = [layer_name]

        # add the layers, using the _add_layer() method
        for i in range(len(layer_path)):
            self._add_layer(
                layer_path[i],
                layer_name[i],
                replace=replace,
                trim=trim,
                src_crs=src_crs,
                dest_crs=dest_crs,
            )

    def remove_layer(self, name: str):
        """Remove a layer in the mask given layer name.

        Args:
            name (str): The name of the layer to be removed.
        """
        if name in self.layers:
            self.layers[name].close()
            del self.layers[name]
        else:
            raise KeyError(f"No layer name {name} found in the mask.")
        self.saved = False

    def rename_layer(self, layer_name: str, new_name: str):
        """Rename a layer in the mask given name.

        Args:
            layer_name (str): The name of the layer to be renamed.
            new_name (str): The new name of the layer.

        Raises:
            KeyError: If the layer name is not found in the mask.
        """

        if layer_name in self.layers:
            self.layers[new_name] = self.layers[layer_name]
            del self.layers[layer_name]
            self.saved = False
        else:
            raise KeyError(f"No layer name {layer_name} found in the mask.")

    def get_res(
        self, layers: Optional[Iterable[str]] = None, product: bool = False
    ) -> dict[str, tuple[float, float]]:
        """
        Get resolutions of layers.

        Args:
            layers (Iterable[str]): The list of layers to get resolution. If not specified,
                select all layers in the object.
            product (bool): Whether the method returns the product of the grid cell height and width.
                By default this is set to False.

        Returns:
            dict: A dictionary with layer names as keys and resolution as values.
        """

        layers = {k: self.layers[k] for k in layers} if layers else self.layers
        res: dict[tuple[float, float]] = {k: val.res for k, val in layers.items()}

        return {k: (res[k][0] * res[k][1]) for k in res} if product else res

    def get_bounds(
        self, layers: Optional[str] = None
    ) -> dict[str, tuple[float, float, float, float]]:
        """Get boundary information of layers

        Args:
            layers (list): The list of layers to get bounds.
                If not specified, select all layers in the object.

        Returns:
            dict: A dictionary with layer names as keys and bounds as values.
        """

        layers = {k: self.layers[k] for k in layers} if layers else self.layers
        return {k: val.bounds for k, val in layers.items()}

    def find_inner_bound(
        self, layers: Optional[str] = None
    ) -> tuple[float, float, float, float]:
        """Find the bounding box that's within specified layers.

        Args:
            layers (list): The list of layers to find inner bound.
                If not specified, select all layers in the object.

        Returns:
            tuple[float, float, float, float]: A bounding box `(minx, miny, maxx, maxy)`
            that represents the bounding box that bound the joint region.
        """

        layer_bound = self.get_bounds(layers=layers)
        return shapely.intersection_all(
            [shapely.box(*bound) for bound in layer_bound.values()]
        ).bounds

    def crop_layer(
        self,
        layer_name: str,
        bounds: tuple[float, float, float, float],
        lat_lon_bounds: bool = True,
        dest_layer_name: Optional[str] = None,
    ):
        """Crop a layer of mask object using `crop_raster` method.

        Args:
            layer_name (str): The name of the layer to be cropped.
            bounds (tuple): The bounds of the cropped layer.
                The order of the bounds should be (left, right, bottom, top).
            lat_lon_bounds (bool): Whether the bounds are in latitude/longitude coordinates.
                By default this is set to True.
            dest_layer_name (str): The name of the cropped layer.
                If not specified, the cropped layer will replace the original layer.
        """

        dest_layer_name = layer_name if dest_layer_name is None else dest_layer_name
        self.layers[dest_layer_name] = crop_raster(
            self.layers[layer_name], bounds, lat_lon_bounds
        )
        self.saved = False

    def trim_layer(self, layer_name: str, dest_layer_name: Optional[str] = None):
        """Trim a layer of mask object using `trim_raster` method.

        Args:
            layer_name (str): The name of the layer to be trimmed.
            dest_layer_name (str): The name of the trimmed layer.
                If not specified, the trimmed layer will replace the original layer.
        """
        dest_layer_name = layer_name if dest_layer_name is None else dest_layer_name
        self.layers[dest_layer_name] = trim_raster(self.layers[layer_name])
        self.saved = False

    def filter_layer(
        self,
        layer_name: str,
        dest_layer_name: Optional[str] = None,
        values: Optional[Iterable[float]] = None,
        min_bound: Optional[float] = None,
        max_bound: Optional[float] = None,
        binarize: bool = False,
    ):
        """Filter a layer of mask object using `filter_raster` method.

        Args:
            layer_name (str): The name of the layer to be filtered.
            dest_layer_name (str): The name of the filtered layer. If not specified,
                the filtered layer will replace the original layer.
            values (list): The list of numberic values in the raster array to filter out as 1.
                If not specified, the method will not filter out any values in the raster.
            min_bound (float): The minimum value in the raster array to filter out as 1.
                If not specified, the method will not filter out any values lower than the min_bound.
            max_bound (float): The maximum value in the raster array to filter out as 1
                If not specified, the method will not filter out any values higher than the max_bound.
            binarize (bool): Whether the method will return 0 and 1 for the raster

        Raises:
            ValueError: If no values, min_bound, nor max_bound are specified.
        """
        if not dest_layer_name:
            dest_layer_name = layer_name

        self.layers[dest_layer_name] = filter_raster(
            self.layers[layer_name], values, min_bound, max_bound, binarize
        )
        self.saved = False

    def merge_layer(
        self,
        method: MergingMethods = "and",
        weights: Optional[dict[str, float]] = None,
        layers: Optional[Iterable[str]] = None,
        trim: bool = False,
        show_raster: bool = True,
        reference_layer: Optional[str] = None,
        attribute_save: bool = True,
        **kwargs,
    ) -> ras.DatasetReader:
        """Merge multiple and flatten multiple layers from self.layers using either 'and'
        and 'sum' method. By default, 'and' method is used, and we would save the result
        to `merged_mask` attribute.

        Adjusted from rasterio's documentation on merging:
        https://rasterio.readthedocs.io/en/latest/api/rasterio.merge.html

        Geospatial bounds and resolution of a new output file in the units of the input
        file coordinate reference system may be provided, but by default, the method will
        use the layer with the best resolution for the output bounds/resolution, unless a
        reference layer is provided.

        Args:
            method (str): The method to merge the layers. By default, 'and' method is used.
                'and' method will perform "AND" operation over the overlapping pixels.
                'sum' method will take the sum of the overlapping pixels.
            weights (dict): The weight of each layer to be merged. By default, all layers
                will have the same weight of 1.
            layers (Iterable[str]): The list of layers to be merged. If not specified, all layers
                in the object will be merged.
            trim (bool): Whether the method will trim the all-empty row/column border
                of the raster. By default, this is set to `False`.
            show_raster (bool): Whether the method will plot the merged raster. By default,
                this is set to `True`.
            reference_layer (str): The name of the layer to be used as the reference layer.
                If not specified, the method will use the layer with the best (highest) resolution.
            attribute_save (bool): Whether the method will save the merged raster to the
                `merged_mask` attribute. By default, this is set to `True`.
            **kwargs: Additional arguments passed to the `rasterio.merge` method.

        Returns:
            rasterio.DatasetReader: The merged raster.

        Raises:
            ValueError: If the specified method is other than 'and' and 'sum'.
        """

        if not layers:
            layers = list(self.layers)

        # specify layers to be used in the merged process
        layers = {k: self.layers[k] for k in layers}

        # If no reference layer is provided, using the layer with best (highest) resolution
        if not reference_layer:
            resolutions = {k: np.prod(v.res) for k, v in layers.items()}
            reference_layer = min(resolutions, key=resolutions.get)

        if method == "and":
            # make sure highest resolution layer is the first input for ras.merge
            merging_layers = [layers[reference_layer]]
            layers.pop(reference_layer)
            merging_layers += list(layers.values())

            arr, aff = merge(merging_layers, method=_and_method, **kwargs)

        elif method == "sum":
            if not weights:
                logger.info(
                    "No weight dictionary provided, all the layers for merging will "
                    "have weights of 1 by default"
                )
                weights = {k: 1 for k in layers}
            else:
                for k in layers.keys():
                    if k not in set(weights.keys()):
                        weights[k] = 1

            temp_layers = {
                name: create_temp_tif(
                    layers[name].read(1) * weight, layers[name].transform
                )
                for name, weight in weights.items()
            }

            # make sure highest resolution layer is the first input for ras.merge
            merging_layers = [temp_layers[reference_layer]]
            temp_layers.pop(reference_layer)
            merging_layers += list(temp_layers.values())

            arr, aff = merge(merging_layers, method=_sum_method, **kwargs)
            for layer in temp_layers.values():
                layer.close()
        else:
            raise ValueError(f"Method {method} is not supported.")

        return_ras = create_temp_tif(arr[0], aff)

        if trim:
            bounds = self.find_inner_bound()
            return_ras = crop_raster(return_ras, bounds)

        if show_raster:
            raster_show(return_ras, title="Merged Mask")

        if attribute_save is True:
            if self.merged_mask:
                logger.info("Overwriting current merged_mask.")
            self.merged_mask = return_ras
            logger.info("Merged Mask saved as attribute 'merged_mask'.")

        return return_ras

    def remove_merge_layer(self):
        """Remove the saved merged mask."""
        self.merged_mask = None

    def add_shape_layer(
        self,
        shapes: dict[str, shapely.Geometry],
        reference_layer: Optional[str] = None,
        resolution: Optional[tuple[float, float]] = None,
        combine_name: Optional[str] = None,
        exclude: bool = False,
        buffer: float = 0.0,
        src_crs: str = "EPSG:4326",
        buffer_crs: str = "EPSG:6933",
        dst_crs: str = "EPSG:4326",
        **kwargs,
    ):
        """Add shapes to the mask layers. This is different from shape extractions,
        as we will simply treat one shp file as a layer, instead of grabbing the merged
        mask within that shape. This method take in a dictionary of shapes, a resolution
        of the result raster with that shape, and add the shape to the mask object. Users
        can also use a reference layer that is present in the mask object to avoid
        manuelly finding resolution.

        Args:
            shapes (dict): A dictionary of key, shape pair.
                Shapes should be a supported geometry type in shapely.
            reference_layer (str): Name to the layer which bounds/resolution is used.
                If not specified, the method will use the layer with the best (highest) resolution.
            resolution (tuple[float, float]): A tuple of `(width_resolution, height_resoution)`.
                If specified with a reference layer, the method ignore the resolution of the referenced layer
                and use the input resolution instead.
            combine_name (str): The name of the combined shape. If specified here, Mask will combine all
                input shapes into one layer.
            exclude (bool): Whether we want to exclude the area specified by the shape or
                the area not specified by shape. By default, this is set to `False`.
            buffer (float): Round buffer distance in km^2 extending out the shapes.
                If input is greater than 0, this method will give approximate representation of all points
                within this given distance of the shapes objects.
            src_crs (str): The source raster's CRS, by default it is 'EPSG:4326' lat lon coordinate system.
            buffer_crs (str): The CRS for the buffer. By default it is 'EPSG:6933' meter coordinate system.
            dst_crs (str): The destination CRS, by default it is 'EPSG:4326' lat lon coordinate system.
            **kwargs: Additional arguments passed to the `rasterio.features.geometry_mask` method.
        """

        if reference_layer:
            shape_transform = self.layers[reference_layer].transform
            if not resolution:
                resolution = self.layers[reference_layer].shape
            else:
                logger.info(
                    "The resolution of the new layers(s) will be %s.", resolution
                )

        else:
            if not resolution:
                raise TypeError(
                    "Please provide the tuple of 'resolution = (width, height)' for the resolution "
                    "of the new shape raster."
                )

        # convert CRS for shapes
        if src_crs != dst_crs:
            for key, shp in shapes.items():
                shapes[key] = shp = convert_shape_crs(
                    shp, src_crs=src_crs, dst_crs=dst_crs
                )

        # add buffer
        if buffer > 0:
            for key, shp in shapes.items():
                if not buffer_crs:
                    raise ValueError(
                        "Please specify a CRS value for projecting shapes for buffer."
                    )

                # convert the shape to a CRS with meter as unit, buffer the shape, and convert it back
                shp = convert_shape_crs(shp, src_crs=dst_crs, dst_crs=buffer_crs)
                shp = shp.buffer(buffer * 1000)
                shp = convert_shape_crs(shp, src_crs=buffer_crs, dst_crs=dst_crs)
                shapes[key] = shp

        if combine_name:
            shapes = {combine_name: shapely.ops.unary_union(list(shapes.values()))}

        for key, shp in shapes.items():
            if not isinstance(shp, shapely.geometry.multipolygon.MultiPolygon):
                # make sure that the shape is a geometry.multipolygon
                shp = shapely.geometry.multipolygon.MultiPolygon([shp])

            shape_bounds = ras.features.bounds(shp)
            if not reference_layer:
                shape_transform = ras.transform.from_bounds(*shape_bounds, *resolution)

            arr = ras.features.geometry_mask(
                shp.geoms if SHAPELY_NEW_API else shp,
                out_shape=resolution,
                transform=shape_transform,
                invert=exclude,
                **kwargs,
            )

            # create raster
            shape_raster = create_temp_tif(arr.astype(np.int8), shape_transform)

            self.layers[key] = shape_raster
            logger.info("Layer %s added to the mask %s.", key, self.name)

    def extract_shapes(
        self,
        shapes: dict[str, shapely.Geometry] | gpd.GeoDataFrame,
        layer: Optional[str] = None,
        combine_shape: bool = False,
        combine_name: Optional[str] = None,
        show_raster: bool = False,
        attribute_save: bool = True,
        src_crs: str = "EPSG:4326",
        dst_crs: str = "EPSG:4326",
        **kwargs,
    ) -> dict[str, ras.DatasetReader]:
        """Extract the shapes on the (by default) merged_mask layer, and save the result
        shape rasters as a dictionary in `shape_mask` attribute.

        Args:
            shapes (dict or GeoDataFrame): A dictionary or GeoPanda's dataframe of shapes
            layer (str): The layer to extract shape from. If unspecified, this method will extract the shapes
                from the `merged_mask` layer.
            combine_shape (bool): Whether combine the shapes as one shape, and only
                one layer will be added as a result. This option is set to False by default.
            combine_name (str): The name of the combined shape if combine_shape is True.
            show_raster (bool): Whether to plot the result shapes. True by default.
            attribute_save (bool): if the program will want to save the shapes to the shape_mask attributes.
                True by default.
            src_crs (str): The source tif CRS, by default it is 'EPSG:4326' lat lon coordinate system.
            dst_crs (str): The destination CRS, by default it is 'EPSG:4326' lat lon coordinate system.
            **kwargs: Additional arguments passed to the `rasterio.mask.mask` method.

        Returns:
            dict: A dictionary of extracted shapes with shape names as keys and shape rasters as values.

        Raises:
            ValueError: If the shapes input is not a dictionary or GeoDataFrame.
            ValueError: If the layer is not specified and the merged_mask is not ready.
        """

        # GeoDataFrame -> Dict
        if hasattr(shapes, "_geometry_column_name"):
            shapes = shapes._geometry_column_name.to_dict()
        elif hasattr(shapes, "_name") and shapes._name == "geometry":
            shapes = shapes.to_dict()

        # Convert CRS for each shape
        if src_crs != dst_crs:
            for k, v in shapes.items():
                shapes[k] = convert_shape_crs(
                    shapes[v], src_crs=src_crs, dst_crs=dst_crs
                )

        if not layer:
            if not self.merged_mask:
                raise ValueError(
                    "No merged_mask ready for the mask object. Please create one using `merge_layer` call"
                    ", or use a layer input by sepcifiying layer=layer_name."
                )
            else:
                layer = self.merged_mask
        else:
            layer = self.layers[layer]

        if not isinstance(shapes, dict):
            raise ValueError("shapes input must be dictionary (key-value pair).")

        if combine_shape:
            if not combine_name:
                raise ValueError("Please specify combined shape name.")
            shapes = {combine_name: shapely.ops.unary_union(shapes.values())}

        return_shape = {}

        # loop through shapes to extract and create temp tif raster
        for key, shp in shapes.items():
            if not isinstance(shp, list):
                shp = [shp]
            arr, aff = rmask.mask(layer, shp, **kwargs)
            raster = create_temp_tif(arr[0], aff)
            return_shape[key] = raster
            if attribute_save:
                if key in self.shape_mask:
                    self.shape_mask[key].close()
                    logger.info(
                        "[Overwritten] Extracted shape %s added to attribute 'shape_mask'.",
                        key,
                    )
                else:
                    logger.info(
                        "Extracted shape %s added to attribute 'shape_mask'.", key
                    )

        if show_raster:
            for name in return_shape:
                raster_show(return_shape[name], title=name)

        if attribute_save:
            self.shape_mask.update(return_shape)
            self.saved = False

        return return_shape

    def remove_shapes(self, names: Iterable[str]):
        """Remove a shape from shape_mask dictionary

        Args:
            names (Iterable[str]): The names of the shapes to be removed.
        """

        for name in names:
            if name not in self.shape_mask.values():
                raise KeyError(f"Shape mask {name} not found in the object.")
            self.shape_mask[name].close()
            del self.shape_mask[name]

    def load_merged_xr(self) -> xr.DataArray:
        """Using xarray to load the merged_layer masks.

        Returns:
            xarray.DataArray: The merged mask xarray.

        Raises:
            ValueError: If the Mask object has not been saved.
        """
        if not self.saved:
            raise ValueError(
                "The Mask object has not been saved, please call save_mask() first."
            )

        merge_path = os.path.join(self.mask_dir, self.name, "merged_mask")
        merge_tif_path = os.path.join(merge_path, "merged_mask.tif")

        logger.info(
            "Please remember to close the merged_mask xarray for further changes of the mask object."
        )
        return ras_to_xarr(merge_tif_path)

    def load_shape_xr(self, names: Optional[str] = None) -> dict[str, xr.DataArray]:
        """Using xarray to load the shape masks.

        Args:
            names (list): The names of the shape masks to be loaded. If not specified, load all shape masks.

        Returns:
            dict: A dictionary of shape masks with shape names as keys and xarray DataArrays as values.

        Raises:
            ValueError: If the Mask object has not been saved.
        """

        if not self.saved:
            raise ValueError(
                "The Mask object has not been saved, please call save_mask() first."
            )

        shape_path = os.path.join(self.mask_dir, self.name, "shape_mask")
        if not names:
            names = self.shape_mask.keys()

        xrs = {
            name: ras_to_xarr(os.path.join(shape_path, f"{name}.tif")) for name in names
        }
        logger.info(
            "Please close the shape_mask xarray(s) for further changes of the mask object."
        )

        return xrs

    def close_files(self):
        """Close all the opened rasters. This method will disable further save_mask() call."""

        for layer in self.layers.values():
            layer.close()

        if self.merged_mask:
            self.merged_mask.close()

        if self.shape_mask:
            for mask in self.shape_mask.values():
                mask.close()

    def save_mask(
        self,
        name: Optional[str] = None,
        mask_dir: Optional[str] = None,
        close_files: bool = False,
    ):
        """Save a mask object to the directory. By default, the directory should be the
        MASK_DIR in config.py.

        **Note**: It is recommended that multiple masks would not be opened
        at the same time to avoid permission issues.

        Args:
            name (str): Name of the mask object to be saved.
                If not specified, use the name stored in this object.
            mask_dir (str): Path to save the mask to. By default, use the value in mask_dir attribute.
            close_files (bool): Whether the program will close all the opened raster files after saving them.
                This will disable further save_mask() call. False by default.
        """

        if self.saved:
            logger.info("Mask %s successfully saved at %s", self.name, self.mask_dir)
            return

        # if user want to use another name to save the mask
        if not name:
            name = self.name
        extension = ".tif"

        if not mask_dir:
            mask_dir = self.mask_dir

        # Ensure mask directory exists
        if not os.path.isdir(mask_dir):
            os.mkdir(mask_dir)

        # get mask object path
        obj_path = os.path.join(mask_dir, self.name)
        if not os.path.isdir(obj_path):
            os.mkdir(obj_path)

        # layer path
        layer_path = os.path.join(obj_path, "layers")
        if not os.path.isdir(layer_path):
            os.mkdir(layer_path)
        else:  # if directory do not exist, remove all the previously saved content
            for f in os.listdir(layer_path):
                if f.split(".")[0] not in self.layers:
                    os.remove(os.path.join(layer_path, f))

        # update layer
        if self.layers:
            for name, raster in self.layers.items():
                save_opened_raster(raster, os.path.join(layer_path, f"{name}.tif"))
                if not close_files:
                    self.layers[name] = ras.open(
                        os.path.join(layer_path, f"{name}.tif")
                    )

        # merge_layer path
        merge_path = os.path.join(obj_path, "merged_mask")
        merge_tif_path = os.path.join(merge_path, "merged_mask.tif")

        if self.merged_mask:
            if not os.path.isdir(merge_path):
                os.mkdir(merge_path)
            save_opened_raster(self.merged_mask, merge_tif_path)
            if not close_files:
                self.merged_mask = ras.open(merge_tif_path)
        else:
            if os.path.isfile(merge_tif_path):  # deleted previously one
                os.remove(merge_tif_path)

        # shape path
        shape_path = os.path.join(obj_path, "shape_mask")
        if not os.path.isdir(shape_path):
            os.mkdir(shape_path)
        else:
            for f in os.listdir(shape_path):
                if f.split(".")[0] not in self.shape_mask:
                    os.remove(os.path.join(shape_path, f))
        if self.shape_mask:
            for name, raster in self.shape_mask.items():
                save_opened_raster(raster, os.path.join(shape_path, name + extension))
                if not close_files:
                    self.shape_mask[name] = ras.open(
                        os.path.join(shape_path, name + extension)
                    )

        self.saved = True
        logger.info("Mask %s successfully saved at %s", self.name, self.mask_dir)

    @classmethod
    def from_name(cls, name: str, mask_dir: str = MASK_DIR):
        """Load a previously saved mask object by name.

        Args:
            name (str): Name of the saved mask.
            mask_dir (str): Directory to look for previously saved mask file and where the mask
                object will be updated. By default, it should be the MASK_PATH specified in config.py

        Returns:
            Mask: The loaded mask object.
        """
        obj_dir = os.path.join(mask_dir, name)

        # create new mask object
        prev_mask = cls(name=os.path.basename(obj_dir))
        prev_mask.mask_dir = mask_dir

        # load the layers
        layer_path = os.path.join(obj_dir, "layers")
        layer_lst = []
        if os.path.isdir(layer_path):
            for i in os.listdir(layer_path):
                prev_mask.layers[i.split(".")[0]] = ras.open(
                    os.path.join(layer_path, i)
                )
                layer_lst.append(i.split(".")[0])
        else:
            logger.info("No previously saved mask found for mask %s", name)

        logger.info("Layer %s loaded to the mask %s.", layer_lst, prev_mask.name)

        # load the merged layer, if there is one
        merge_path = os.path.join(obj_dir, "merged_mask")
        merge_file_path = os.path.join(merge_path, "merged_mask.tif")
        if os.path.exists(merge_file_path):
            prev_mask.merged_mask = ras.open(merge_file_path)
            logger.info("Merged_mask loaded to the mask %s.", prev_mask.name)
        else:
            logger.info("No merged Mask found.")

        # load the shape masks
        shape_path = os.path.join(obj_dir, "shape_mask")
        shape_lst = []
        if os.path.isdir(shape_path):
            for i in os.listdir(shape_path):
                prev_mask.shape_mask[i.split(".")[0]] = ras.open(
                    os.path.join(shape_path, i)
                )
                shape_lst.append(i.split(".")[0])
            logger.info(
                "Shape mask %s loaded to the mask %s.", shape_lst, prev_mask.name
            )

        prev_mask.saved = True
        return prev_mask


def open_tif(
    layer_path: str, show_raster: bool = True, close: bool = False
) -> ras.DatasetReader:
    """Open a layer using rasterio given a file path.

    Args:
        layer_path (str): The path to the TIFF file.
        show_raster (bool): Whether the method will plot the raster. True by default.
        close (bool): Whether the method will close the raster afterward, prevent possible permission errors.
            False by default.

    Returns:
        rasterio.DatasetReader: The opened raster.

    Raises:
        FileNotFoundError: If the layer_path does not exist.
    """

    if not os.path.exists(layer_path):
        raise FileNotFoundError(f"{layer_path} not found.")
    openning_tif: ras.DatasetReader = ras.open(layer_path)

    if show_raster:
        raster_show(openning_tif)
    if not close:
        logger.info("Please remember to close the file with .close()")
    else:
        openning_tif.close()

    return openning_tif


def ras_to_xarr(
    raster: str, band_name: Optional[str] = None, adjust_coord: bool = True
) -> xr.DataArray:
    """Open a raster (rasterio openner) with rioxarray

    Args:
        raster (str): The path to the raster
        band_name (str): The name of the band. This is optional. If no band name is specified,
            no band will be renamed.
        adjust_coord (bool): Whether to adjust the coordinate names to `lat`/`lon`. True by default.

    Returns:
        xarray.DataArray: The opened raster as xarray DataArray.
    """

    xarr: xr.DataArray = open_rasterio(raster)

    if adjust_coord:
        xarr = xarr.rename({"x": "lon", "y": "lat"})
        xarr = xarr.sortby(["lat", "lon"])
    if band_name:
        xarr = xarr.rename({"band": band_name})

    return xarr


def create_temp_tif(
    arr: np.ndarray, transform: ras.Affine, open_raster: bool = True
) -> ras.DatasetReader | str:
    """Create a ras.DatasetReader object openning a temporary rasterio file

    **Note**: Arguments `arr` and `transform` input can be the output
    for `ras.merge.merge` and `ras.mask.mask`.

    Args:
        arr (ArrayLike): An ArrayLike object that contains values of the layer.
        transform (rasterio.Affine): Affine transformation for the layer.
        open_raster (bool): Whether the raster will be opened. True by default.

    Returns:
        rasterio.DatasetReader: The temporary raster.
    """

    with MemoryFile() as memfile:
        with ras.open(
            memfile.name,
            "w",
            driver="GTiff",
            height=arr.shape[0],
            width=arr.shape[1],
            count=1,
            dtype=arr.dtype,
            compress="lzw",
            crs="+proj=latlong",
            transform=transform,
        ) as dst:
            dst.write(arr, 1)

        if open_raster:
            return ras.open(memfile.name)

        return memfile.name


def save_opened_raster(raster: ras.DatasetReader, path: str):
    """Helper method to close the opened raster and save it to the given path.

    Args:
        raster (rasterio.DatasetReader): The opened raster.
        path (str): The path to save the raster.
    """

    arr, transform = raster.read(1), raster.transform
    raster.close()
    save_raster(arr, transform, path)


def save_raster(arr: np.ndarray, transform: ras.Affine, path: str):
    """Given an array and an Affine transform, save the raster as tif file to the specified path.
    https://rasterio.readthedocs.io/en/latest/topics/writing.html

    **Note**: Arguments `arr` and `transform` input can be the output
    for `ras.merge.merge` and `ras.mask.mask`.

    Args:
        arr (ArrayLike): An array with values of the raster.
        transform (rasterio.Affine): Affine transform information for that layer.
        path (str): the path to where the tif file is saved.
    """

    with ras.open(
        path,
        "w",
        driver="GTiff",
        height=arr.shape[0],
        width=arr.shape[1],
        count=1,
        dtype=arr.dtype,
        compress="lzw",
        crs="+proj=latlong",
        transform=transform,
    ) as dst:
        dst.write(arr, 1)


def crop_raster(
    raster: ras.DatasetReader,
    bounds: tuple[float, float, float, float],
    lat_lon_bounds: bool = True,
) -> ras.DatasetReader:
    """Crop a raster given geographical coordinates or its array index bounds

    Args:
        raster (ras.DatasetReader): An opened rasterio tif dataset.
        bounds (tuple[float, float, float, float]): Geographical coordinates or array index bounds.
        lat_lon_bounds (bool): If True (by default), the bounds tuple are latitude-longitude coords input.

    Returns:
        ras.DatasetReader: The cropped raster.
    """

    if lat_lon_bounds:  # if bounds are lat/long
        window = ras.windows.from_bounds(
            bounds[0], bounds[1], bounds[2], bounds[3], raster.transform
        )
    else:  # bounds are index
        window = ras.windows.Window.from_slices(
            (bounds[0], bounds[1]), (bounds[2], bounds[3])
        )

    with MemoryFile() as memfile:
        kwargs = raster.meta.copy()
        kwargs.update(
            {
                "height": window.height,
                "width": window.width,
                "transform": ras.windows.transform(window, raster.transform),
            }
        )

        with ras.open(memfile.name, "w", compress="lzw", **kwargs) as dst:
            dst.write(raster.read(window=window))
            dst.close()

        return ras.open(memfile.name)


def reproject_raster(
    src: ras.DatasetReader,
    src_crs: str,
    dst_crs: str = "EPSG:4326",
    trim: bool = False,
    **kwargs,
) -> ras.DatasetReader:
    """Convert CRS of TIFF file from one to another. By default, we want to make the destination raster
    in lat/lon coordinate system. However, this can still be changed by specifying the dst_crs.

    Args:
        src (ras.DatasetReader): The source raster
        src_crs (str): The source CRS
        dst_crs (str): By default, we want to make the destination raster in lat/lon coordinate system.
            However, this can still be changed by specifying the dst_crs.
        trim (bool): False by default. If True, we will trim the empty border to save space.

    Returns:
        ras.DatasetReader: The reprojected raster.
    """

    transform, width, height = ras.warp.calculate_default_transform(
        src_crs, dst_crs, src.width, src.height, *src.bounds
    )
    kwargs = src.meta.copy()
    kwargs.update(
        {"crs": dst_crs, "transform": transform, "width": width, "height": height}
    )

    # write it to another file: the CRS corrected one
    # rasterio.readthedocs.io/en/latest/topics/reproject.html
    with MemoryFile() as memfile:
        with ras.open(memfile.name, "w", compress="lzw", **kwargs) as dst:
            for i in range(1, src.count + 1):
                ras.warp.reproject(
                    source=ras.band(src, i),
                    destination=ras.band(dst, i),
                    src_transform=src.transform,
                    src_crs=src_crs,
                    dst_transform=transform,
                    dst_crs=dst_crs,
                    resampling=ras.warp.Resampling.nearest,
                )

        logger.info("Raster %s has been reprojected to %s CRS.", src.name, dst_crs)
        return_ras = ras.open(memfile.name)

        if trim:
            return trim_raster(return_ras)

        return return_ras


def apply_fn_to_raster(raster: ras.DatasetReader, fn: callable):
    """Apply a custom function to individual values in the raster array.

    Args:
        raster (ras.DatasetReader): The source raster.
        fn (callable): A function/method to be applied to the raster array
            the function needs map individual values to another.
            This callable will be JIT compiled, if possible.
        kwargs: the arguments to be passed to the function

    return: (ras.DatasetReader) The updated raster
    """

    orig_values = raster.read(1)
    fn = njit(fn)

    @njit(parallel=True)
    def _dummy_vectorize_fn(fn: callable, arr: np.ndarray):
        assert arr.ndim == 2

        for i in prange(arr.shape[0]):
            for j in prange(arr.shape[1]):
                arr[i, j] = fn(arr[i, j])

        return arr

    return create_temp_tif(_dummy_vectorize_fn(fn, orig_values), raster.transform)


def filter_raster(
    raster: ras.DatasetReader,
    values: Optional[Iterable[float]] = None,
    min_bound: Optional[float] = None,
    max_bound: Optional[float] = None,
    binarize: bool = False,
) -> ras.DatasetReader:
    """Filter the raster with a list of values to be True. This method can also set
    a min bound and max bound to select values in the raster higher or lower than
    these two boundaries. If binarize is False (by default), the method will return
    the original values of the raster that satisfy the condition.

    Args:
        raster (ras.DatasetReader): The source raster
        values (list): The list of numberic values in the raster array to filter out as 1.
            Optional. If not specified, the method will not filter out any values in the raster.
        min_bound (float): The minimum value in the raster array to filter out as 1.
            Optional. If not specified, the method will not filter out any values lower than the min_bound.
        max_bound (float): The maximum value in the raster array to filter out as 1
            Optional. If not specified, the method will not filter out any values higher than the max_bound.
        binarize (bool): Whether the method will return 0 and 1 for the raster

    Returns:
        ras.DatasetReader: The filtered raster.

    Raises:
        ValueError: If no values, min_bound, nor max_bound are specified.
    """
    if values is None and min_bound is None and max_bound is None:
        raise ValueError("Please specify any of values, min_bound, or max_bound.")

    bool_arr = raster.read(1)
    if values is not None:
        bool_arr = np.isin(bool_arr, values)
    if min_bound is not None:
        bool_arr = bool_arr > min_bound
    if max_bound is not None:
        bool_arr = bool_arr < max_bound

    # if the method return 0 and 1 for the raster
    if binarize:
        return create_temp_tif(bool_arr.astype(np.uint8), raster.transform)
    return create_temp_tif(bool_arr * raster.read(1), raster.transform)


def trim_raster(raster: ras.DatasetReader) -> ras.DatasetReader:
    """Remove the all-zero columns and rows at the border of the raster and returns the trimmed raster.
    This method does not remove the all-zero cols or rows in the middle of the valid values.

    For example: if the raster.read(1) (the array values) is the following::

        np.array(
            [
                [0, 0, 9, 0, 0, 0, 0],
                [0, 0, 1, 2, 0, 0, 0],
                [0, 0, 2, 3, 4, 0, 0],
                [0, 0, 0, 0, 0, 0, 0],
                [0, 0, 1, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0, 0],
            ]
        )

    The bounds of valid values will be (2, 4, 0, 4);
    The all zero columns and rows at the border of the array will be removed.

    Args:
        raster (ras.DatasetReader): The source raster

    Returns:
        ras.DatasetReader: The trimmed raster.
    """

    arr = raster.read(1)

    all_zero_col = np.argwhere(np.all(arr[..., :] == 0, axis=0)).reshape(-1)
    all_zero_row = np.argwhere(np.all(arr[:, ...] == 0, axis=1)).reshape(-1)

    l_idx, t_idx = 0, 0
    r_idx, b_idx = arr.shape[1] - 1, arr.shape[0] - 1

    # find left and right index where all empty column start and ends
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

    # find top and bottom index where all empty column start and ends
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

    return crop_raster(raster, bounds, lat_lon_bounds=False)


def filter_area(
    self,
    min_area: float,
    layer_name: Optional[str] = None,
    dest_layer_name: Optional[str] = None,
    shape_value: int = 1,
    src_crs: str = "EPSG:4326",
    area_calc_crs: str = "EPSG:6933",
    connectivity: Literal[4, 8] = 4,
) -> ras.DatasetReader:
    """Eliminate the small area of a certain value in the raster by converting it to
    an equal-area reprojected series of connected shapes, removing shapes that are
    smaller than the given area in km^2, and turning it back to a raster.

    Args:
        min_area (float): The area threshold in km^2. Any connected group with smaller area
            then this parameter will be removed from the raster.
        layer_name (str): The name of the raster to be selected from the mask object
            If not specified, the method will use the `merged_mask` layer
        dest_layer_name (str): The name of the new raster to be formed from eliminating
            small area in the original layer. If not specified, this function will
            return the raster instead of assigning it to the layers.
        shape_value (int): The value of the grid cells to be groups for elimination.
            1 by default. (Find all the connected groups of cells with area 1)
        src_crs (str): The CRS of the raster to be selected from the mask object.
            EPSG:4326 by default.
        area_calc_crs (str): The CRS used for reprojecting the raster for area calculation.
            EPSG:6933 by default. (Equal area projection)
        connectivity (int): Should either be 4 or 8. If 4, then the 4 surrounding cells of
            each grid cell will be counted for grouping connected grid cells, otherwise, all
            8 surround cells of each grid cell will be counted.

    Returns:
        ras.DatasetReader: The raster with small area eliminated.

    Raises:
        KeyError: If the layer_name does not exist in the mask object.
    """

    if layer_name is None:
        raster = self.merged_mask
    elif layer_name in self.layers:
        raster = self.layers[layer_name]
    else:
        raise KeyError("Specified raster does not exist in the mask object.")

    logger.info(
        "Reprojecting the raster to %s for equal area calculation.", area_calc_crs
    )
    # reproject the raster into equal area projection
    reproj = reproject_raster(raster, src_crs=src_crs, dst_crs=area_calc_crs)

    # general shapes
    res_shapes = ras.features.shapes(
        reproj.read(1).astype("uint8"),
        connectivity=connectivity,
        transform=reproj.transform,
    )

    res_shapes = [
        shapely.geometry.shape(i[0]) for i in res_shapes if i[1] == shape_value
    ]

    g_series = gpd.GeoSeries(res_shapes)

    # filter the shapes by area
    filtered_shape = g_series[g_series.area > min_area * (1000 * 1000)]

    # set up transform and shape for shape to raster
    shape_transform = raster.transform
    resolution = raster.shape

    shp = shapely.geometry.multipolygon.MultiPolygon(filtered_shape.values)
    shp = convert_shape_crs(shp, src_crs=area_calc_crs, dst_crs=src_crs)

    logger.info("Reverting the remaining shapes back to a raster.")
    arr = ras.features.geometry_mask(
        shp.geoms if SHAPELY_NEW_API else shp,
        out_shape=resolution,
        transform=shape_transform,
        invert=True,
    )

    # create raster
    shape_raster = create_temp_tif(arr.astype(np.int8), shape_transform)

    if dest_layer_name:
        self.layers[dest_layer_name] = shape_raster

    return shape_raster


def convert_shape_crs(shape: shapely.Geometry, src_crs: str, dst_crs: str):
    """Convert the CRS of a shape object given its source CRS and destinated CRS.

    Args:
        shape (shapely.Geometry): The shape object to be converted.
        src_crs (str): The source CRS of the shape object.
        dst_crs (str): The destinated CRS of the shape object.

    Returns:
        shapely.Geometry: The converted shape object.
    """

    src_crs = pyproj.CRS(src_crs)
    dst_crs = pyproj.CRS(dst_crs)

    buffer_project = pyproj.Transformer.from_crs(
        src_crs, dst_crs, always_xy=True
    ).transform
    return shapely.ops.transform(buffer_project, shape)


def _sum_method(merged_data, new_data, merged_mask, new_mask, index, roff, coff):
    """The sum method will add up the values from all the layers. We can also
    customize the weights. The behind scene of this method is that it multiplys
    each layers with the corresponding weight, and add the in-memory temporary
    layers together."""

    if len(np.unique(merged_data)) == 1:
        mask = np.ones(merged_data.shape, dtype=bool)
        np.add(
            np.zeros(merged_data.shape),
            new_data.data,
            out=merged_data,
            where=mask,
            casting="unsafe",
        )
    else:
        np.add(merged_data, new_data.data, out=merged_data, casting="unsafe")


def _and_method(merged_data, new_data, merged_mask, new_mask, index, roff, coff):
    """By default, the merge_layer method will use a binary 'and' method:
    if any of the n grid cells of the n layers at the same location have 0,
    then the returned self.merged_layer will also have 0 at that location.
    In other words, if all the layers indicate that a land is not unavailable (!=0),
    the merged result will have value 1."""

    if len(np.unique(merged_data)) == 1:
        mask = np.ones(merged_data.shape, dtype=bool)
    else:
        mask = merged_data != 0

    np.copyto(merged_data, new_data.data, where=mask, casting="unsafe")


def show(
    raster: ras.DatasetReader | str,
    shape: Optional[gpd.GeoSeries] = None,
    shape_width: float = 0.5,
    shape_color: str = "black",
    figsize: tuple[float, float] = (10, 6),
    title: Optional[str] = None,
    title_size: float = 12,
    lat_lon: bool = True,
    colorbar: bool = True,
    grid: bool = False,
    return_fig: bool = False,
    **kwargs,
) -> Optional[plt.Figure]:
    """Plot a raster file.

    Args:
        raster (ras.DatasetReader | str): A `rasterio.DatasetReader` or string. If a string is passed,
            the method will open the raster file using `open_tif()` method.
        shape (geopandas.GeoSeries): Shapes to be plotted over the raster. None by default.
        shape_width (float): The line width for plotting shapes. 0.5 by default.
        shape_color (str): Color of the shape border. Black by default.
        figsize (tuple[float, float]): The size of the figure. (10, 6) by default.
        title (str): The title of the plot. If not specified, the title will be the name of the raster.
        lat_lon (bool): Whether the program will show the lat/long coordinates in the plot. True by default.
        colorbar (bool): Whether the program shows the legend for the values of the plot. True by default.
        grid (bool): whether the program shows the grid. False by default.
        return_fig (bool): Whether to return the matplotlib figure. False by default.
        **kwargs: Other keyword arguments to be passed to `matplotlib.pyplot.imshow()`.

    Returns:
        matplotlib.figure.Figure: The matplotlib figure. Only returned if `return_fig` is True.
    """

    if isinstance(raster, str):
        raster = open_tif(raster, show_raster=False)

    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot()

    if shape is not None:
        shape.boundary.plot(ax=ax, linewidth=shape_width, color=shape_color)

    if lat_lon:
        ax_img = ax.imshow(
            raster.read(1),
            interpolation="none",
            extent=plotting_extent(raster.read(1), raster.transform),
            **kwargs,
        )
    else:
        ax_img = ax.imshow(raster.read(1), interpolation="none", **kwargs)

    if title is True:
        title = raster.name
    ax.set_title(title, size=title_size)
    if colorbar:
        uniques = np.unique(raster.read(1))
        if len(uniques) < 3:
            plt.colorbar(ax_img, ax=ax, ticks=uniques)
        else:
            plt.colorbar(ax_img, ax=ax)
    if grid:
        ax.grid()
    if return_fig:
        return fig


raster_show = show


def show_all(r_dict: dict[str, ras.DatasetReader], **kwargs):
    """Show all the layers given a dictionary of `name: raster` pairs.
    Example use case: show_all(my_mask.layers); show_all(my_mask.shape_mask).

    Args:
        r_dict (dict): A dictionary of `name: raster` pairs.
        **kwargs: Other keyword arguments to be passed to `show()` method.
    """
    for name, r in r_dict.items():
        show(r, title=name, **kwargs)
