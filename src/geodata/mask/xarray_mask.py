"""Xarray-native mask adapter for applying saved Mask objects to datasets."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Literal, cast

import numpy as np
import xarray as xr

from .spatial import calc_grid_area, coarsen, ds_reformat_index


def _ensure_dataset(data: xr.Dataset | xr.DataArray) -> xr.Dataset:
    if isinstance(data, xr.Dataset):
        return data
    name = data.name or "value"
    return data.to_dataset(name=name)


def _to_mask_2d(mask: xr.DataArray) -> xr.DataArray:
    mask = mask.reset_coords(drop=True)
    if "band" in mask.dims:
        mask = mask.isel(band=0, drop=True)
    return mask.transpose("lat", "lon")


@dataclass
class XarrayMask:
    """Mask adapter that aligns saved mask rasters to a target xarray grid."""

    grid: xr.Dataset
    merged_mask: xr.DataArray | None = None
    shape_masks: dict[str, xr.DataArray] = field(default_factory=dict)

    @classmethod
    def from_mask(
        cls,
        mask,
        grid: xr.Dataset | xr.DataArray,
        include_merged: bool = True,
        include_shapes: bool = True,
    ) -> "XarrayMask":
        grid_ds = ds_reformat_index(_ensure_dataset(grid))
        grid_ds = cast(xr.Dataset, grid_ds)
        merged = None
        shapes: dict[str, xr.DataArray] = {}

        if include_merged and mask.merged_mask:
            merged = coarsen(mask.load_merged_xr(), grid_ds)
        if include_shapes and mask.shape_mask:
            shapes = {k: coarsen(v, grid_ds) for k, v in mask.load_shape_xr().items()}

        if merged is None and not shapes:
            raise ValueError(
                f"No mask found in {mask.name}. Please create a proper mask object first."
            )

        return cls(grid=grid_ds, merged_mask=merged, shape_masks=shapes)

    @classmethod
    def from_name(
        cls,
        name: str,
        grid: xr.Dataset | xr.DataArray,
        mask_dir: str | None = None,
        include_merged: bool = True,
        include_shapes: bool = True,
    ) -> "XarrayMask":
        from geodata import Mask  # lazy import to avoid circular imports

        if mask_dir is None:
            from geodata import config

            mask = Mask.from_name(name, mask_dir=config.MASK_DIR)
        else:
            mask = Mask.from_name(name, mask_dir=mask_dir)
        return cls.from_mask(
            mask,
            grid=grid,
            include_merged=include_merged,
            include_shapes=include_shapes,
        )

    @staticmethod
    def compute_grid_area(grid: xr.Dataset | xr.DataArray) -> xr.DataArray:
        xr_ds = ds_reformat_index(_ensure_dataset(grid))
        area_arr = np.zeros((xr_ds.lat.shape[0], xr_ds.lon.shape[0]))
        lat_diff = np.abs((xr_ds.lat[1].values - xr_ds.lat[0].values))
        for i, lat in enumerate(xr_ds.lat.values):
            lat_bottom = lat - lat_diff / 2
            lat_top = lat + lat_diff / 2
            area_arr[i] = np.round(
                calc_grid_area(
                    [
                        (xr_ds.lon.values[0], lat_top),
                        (xr_ds.lon.values[0], lat_bottom),
                        (xr_ds.lon.values[1], lat_bottom),
                        (xr_ds.lon.values[1], lat_top),
                    ]
                ),
                2,
            )
        return xr.DataArray(
            area_arr,
            dims=("lat", "lon"),
            coords={"lat": xr_ds.lat.values, "lon": xr_ds.lon.values},
            name="area",
        )

    def _target_masks(self) -> dict[str, xr.DataArray]:
        res: dict[str, xr.DataArray] = {}
        if self.merged_mask is not None:
            res["merged_mask"] = _to_mask_2d(self.merged_mask)
        for key, value in self.shape_masks.items():
            res[key] = _to_mask_2d(value)
        return res

    def attach(
        self, dataset: xr.Dataset | xr.DataArray, include_area: bool = True
    ) -> dict[str, xr.Dataset]:
        ds = ds_reformat_index(_ensure_dataset(dataset))
        if "time" in ds.dims:
            ds = ds.transpose("time", "lat", "lon")

        masks = self._target_masks()
        if not masks:
            raise ValueError("No masks available in XarrayMask.")

        area = self.compute_grid_area(self.grid) if include_area else None
        out: dict[str, xr.Dataset] = {}
        for key, mask in masks.items():
            cur = ds.assign({"mask": mask})
            if area is not None:
                cur = cur.assign({"area": area})
            out[key] = cur
        return out

    def apply(
        self,
        dataset: xr.Dataset | xr.DataArray,
        mode: Literal["where", "multiply"] = "where",
        include_area: bool = False,
    ) -> dict[str, xr.Dataset]:
        ds = ds_reformat_index(_ensure_dataset(dataset))
        if "time" in ds.dims:
            ds = ds.transpose("time", "lat", "lon")

        masks = self._target_masks()
        if not masks:
            raise ValueError("No masks available in XarrayMask.")
        if mode not in {"where", "multiply"}:
            raise ValueError("mode can only be 'where' or 'multiply'")

        area = self.compute_grid_area(self.grid) if include_area else None
        out: dict[str, xr.Dataset] = {}
        for key, mask in masks.items():
            valid = mask > 0
            cur = ds.copy()
            for var in list(cur.data_vars):
                da = cur[var]
                if "lat" in da.dims and "lon" in da.dims:
                    if mode == "where":
                        cur = cur.assign({var: cast(Any, da.where(valid))})
                    else:
                        cur = cur.assign({var: cast(Any, da * valid)})
            if include_area and area is not None:
                cur = cur.assign({"area": area})
            out[key] = cast(xr.Dataset, cur)
        return out


__all__ = ["XarrayMask"]
