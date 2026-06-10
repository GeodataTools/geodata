"""Mask package namespace.

This package hosts mask-related modules (e.g. spatial helper utilities) while
preserving backward-compatible access to the legacy ``geodata.mask`` module API.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path

from .spatial import calc_grid_area, calc_shp_area, coarsen, ds_reformat_index
from .xarray_mask import XarrayMask

_LEGACY_MODULE_PATH = Path(__file__).resolve().parent.parent / "mask.py"
_LEGACY_SPEC = importlib.util.spec_from_file_location(
    "geodata._legacy_mask_module", _LEGACY_MODULE_PATH
)
if _LEGACY_SPEC is None or _LEGACY_SPEC.loader is None:
    raise ImportError(f"Could not load legacy mask module from {_LEGACY_MODULE_PATH}")
_legacy_mask_module = importlib.util.module_from_spec(_LEGACY_SPEC)
_LEGACY_SPEC.loader.exec_module(_legacy_mask_module)

# Re-export all public names from legacy ``mask.py``.
for _name in dir(_legacy_mask_module):
    if _name.startswith("_"):
        continue
    globals()[_name] = getattr(_legacy_mask_module, _name)

# Keep explicit access to phase-1 extracted helpers in this namespace.
globals().update(
    {
        "ds_reformat_index": ds_reformat_index,
        "coarsen": coarsen,
        "calc_grid_area": calc_grid_area,
        "calc_shp_area": calc_shp_area,
        "XarrayMask": XarrayMask,
    }
)

__all__ = [name for name in globals() if not name.startswith("_")]
