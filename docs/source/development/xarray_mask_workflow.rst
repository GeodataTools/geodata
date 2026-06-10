Xarray masking workflow
=========================

.. note::

   **Audience:** contributors and maintainers. This page records the xarray-first
   masking implementation phases. For usage, see
   :doc:`/mask/xarray_mask_tutorial`.

This page summarizes the **xarray-first masking** work added alongside the
longer-term plan in :doc:`mask_xarray_migration_plan`. The legacy path based on
``Cutout`` (``add_mask``, ``add_grid_area``, ``mask``) is unchanged for now; the
new pieces let you mask **any** model or analysis output
given as an ``xarray.Dataset`` or ``xarray.DataArray``, without threading mask
logic through model classes.

What was added
--------------

**Phase 0 — behavior freeze (tests only)**

Offline tests lock in legacy masking behavior so refactors do not silently change
results:

* Coarsening / alignment of saved mask rasters onto a target grid.
* Grid cell area computation consistent with the cutout-style workflow.
* The structure of outputs from ``Cutout.mask(...)`` (keys, variables, dimensions).
* Selected error paths (missing mask, missing area, invalid mask state).

**Phase 1 — shared spatial helpers**

The following helpers now live in ``geodata.mask.spatial`` and are re-used from
``cutout`` (and plotting code where relevant):

* ``ds_reformat_index`` — normalize coordinates toward ``lat`` / ``lon``.
* ``coarsen`` — align a higher-resolution mask grid to a target grid.
* ``calc_grid_area`` / ``calc_shp_area`` — area utilities used by the masking workflow.

Public names on ``geodata.cutout`` (e.g. ``coarsen``, ``calc_grid_area``) remain
available as aliases for backward compatibility.

**Phase 2 — ``XarrayMask``**

``XarrayMask`` (``from geodata import XarrayMask``) provides:

* ``from_name`` / ``from_mask`` — load a saved ``Mask`` and align
  merged and shape masks to a target ``grid`` (your model output or any dataset
  with compatible ``x``/``y`` or ``lat``/``lon`` coordinates).
* ``compute_grid_area`` — per-cell area on the target grid (same idea as cutout
  grid area).
* ``attach`` — return a dict of datasets like legacy ``Cutout.mask``: original
  variables plus ``mask`` and optional ``area``.
* ``apply`` — return masked data (``mode="where"`` for NaN outside mask,
  ``mode="multiply"`` for zero outside mask), optionally with ``area``.

**Integration pattern (no coupling inside models)**

Masking is intentionally **not** built into wind, pvlib, or other model ``estimate``
APIs. The intended usage is:

1. Run the model and obtain ``output_ds`` (or a ``DataArray`` you wrap in a
   one-variable dataset).
2. Build ``XarrayMask.from_name("my_mask", grid=output_ds, mask_dir=...)`` if needed.
3. Call ``attach(output_ds)`` or ``apply(output_ds, ...)`` for analysis.

See :doc:`/mask/xarray_mask_tutorial` for a step-by-step notebook, and the offline
tests under ``tests/pr/`` (e.g. ``test_xarray_mask.py``,
``test_wind_xarraymask_integration.py``) for concrete examples.

Package layout note
-------------------

The repository currently has both:

* ``src/geodata/mask.py`` — original ``geodata.mask`` implementation (``Mask``,
  raster helpers, etc.).
* ``src/geodata/mask/`` — package namespace that re-exports that API **and**
  hosts new modules (``spatial.py``, ``xarray_mask.py``).

Imports like ``from geodata import Mask`` and ``from geodata import XarrayMask``
continue to work during this transition.

See also
--------

* :doc:`/mask/xarray_mask_tutorial` — step-by-step notebook (offline runnable).
* :doc:`mask_xarray_migration_plan` — full migration phases and deprecation plan.
* :doc:`/legacy/mask_on_cutout` — legacy notebook: masks via ``Cutout``.
* :doc:`/mask/mask_creation_workflow` — building and saving ``Mask`` objects from rasters.
