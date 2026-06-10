# Mask-Without-Cutout Migration Plan

```{note}
**Audience:** contributors and maintainers. For applying saved masks to model
output, see [XarrayMask tutorial](../mask/xarray_mask_tutorial.ipynb).
```

## Goal

Replace Cutout-dependent masking with a direct xarray-based workflow:

- `datasets -> models -> masking -> analysis`

The new masking flow should work on model output (`xarray.Dataset` / `xarray.DataArray`) directly, while reusing current `mask.py` code as much as possible.

## What Changes, What Stays

- Keep:
  - `Mask` object for raster/shapefile mask creation and persistence.
  - Existing layer operations in `src/geodata/mask.py` (`add_layer`, `filter_layer`, `merge_layer`, `extract_shapes`, `save_mask`, `from_name`).
  - Existing geospatial utilities (`ras_to_xarr`, `calc_grid_area` logic from `cutout.py`, coordinate formatting helpers).
- Remove dependency on:
  - `Cutout.add_mask(...)`
  - `Cutout.add_grid_area(...)`
  - `Cutout.mask(...)`
- Add:
  - A new xarray-focused masking adapter class/module (proposed below).

## Proposed Target API

Create a dedicated class (example name: `XarrayMask`) that only deals with xarray data:

1. **Creation / loading**
   - `XarrayMask.from_mask(mask: Mask, grid: xr.Dataset | xr.DataArray, include_merged=True, include_shapes=True)`
   - `XarrayMask.from_name(name: str, grid: xr.Dataset | xr.DataArray, mask_dir=...)`
2. **Area calculation**
   - `XarrayMask.compute_grid_area(grid: xr.Dataset | xr.DataArray) -> xr.DataArray`
3. **Applying mask**
   - `XarrayMask.attach(dataset, include_area=True) -> dict[str, xr.Dataset]`
     - Equivalent to current `Cutout.mask(...)` behavior (mask as extra variables).
   - `XarrayMask.apply(dataset, mode="where", include_area=False) -> dict[str, xr.Dataset]`
     - New convenience method returning mask-applied outputs:
       - `mode="where"`: outside mask -> NaN
       - `mode="multiply"`: outside mask -> 0

This gives both:
- transparent feature-style behavior (`attach`)
- direct filtered outputs (`apply`)

## Reuse Map (Do Not Reinvent)

Directly reuse existing code paths:

- From `src/geodata/mask.py`:
  - `Mask.from_name(...)`
  - `Mask.load_merged_xr()` / `Mask.load_shape_xr()`
- From `src/geodata/cutout.py`:
  - `ds_reformat_index(...)` (move/shared helper)
  - `coarsen(...)` (move/shared helper)
  - `calc_grid_area(...)` (move/shared helper)
- Keep the same coordinate conventions:
  - normalize to `lat`, `lon`
  - align mask grid to target dataset grid before applying

Refactor suggestion:
- Move shared helpers into a new utility module, e.g. `src/geodata/spatial.py` or `src/geodata/mask_xarray.py`, then import from both old and new flows during transition.

## Migration Phases

### Phase 0 - Freeze Current Behavior

- Add tests that lock existing behavior for:
  - coarsening/alignment from mask raster to target grid
  - area computation
  - output structure currently returned by `Cutout.mask(...)`

This prevents regressions while extracting logic.

### Phase 1 - Extract Shared Spatial Helpers

- Move (or duplicate temporarily) these functions out of `cutout.py`:
  - `ds_reformat_index`
  - `coarsen`
  - `calc_grid_area`
- Add unit tests for each helper independent of `Cutout`.

### Phase 2 - Introduce `XarrayMask`

- Implement class that:
  - loads saved `Mask` by name
  - converts mask rasters to xarray
  - coarsens/aligned to target grid
  - computes area from target grid
  - provides `attach()` and `apply()`

### Phase 3 - Integrate into datasets -> models workflow

- At model output point (where xarray result exists), call:
  - `xmask = XarrayMask.from_name("my_mask", grid=model_ds)`
  - `masked = xmask.apply(model_ds, mode="where")`
- Keep `attach()` available for advanced users needing raw mask + area features.

### Phase 4 - Deprecate Cutout Masking Surface

- Mark these as deprecated:
  - `Cutout.add_mask`
  - `Cutout.add_grid_area`
  - `Cutout.mask`
- Keep them as wrappers calling new `XarrayMask` for 1-2 releases.

### Phase 5 - Remove Cutout Dependency

- Remove or archive old mask-coupled Cutout paths once internal usage is migrated.
- Keep `Cutout` only if still needed for data preparation.

## Detailed Behavior Decisions

To avoid ambiguity, define these explicitly:

- Mask value semantics:
  - `mask > 0` means valid/included
  - `mask <= 0` means excluded
- Apply scope:
  - apply to all data variables by default
  - optional include/exclude variable list
- Output keys:
  - `"merged_mask"` for merged mask
  - shape names for shape masks (same as current behavior)
- Alignment:
  - always reformat coords to `lat`/`lon`
  - always transpose to `time, lat, lon` when `time` exists
- Area:
  - computed from target grid only (not from mask grid) to stay consistent with model outputs

## Risks and Mitigations

- Risk: hidden coordinate mismatches (`x/y` vs `lat/lon`, descending latitude).
  - Mitigation: centralize coordinate normalization in one helper and test with both styles.
- Risk: users depending on old `Cutout.mask` output shape.
  - Mitigation: make `attach()` output identical structure and keep temporary wrappers.
- Risk: performance hit when repeatedly coarsening same mask.
  - Mitigation: cache aligned masks keyed by grid signature (lat/lon hashes + mask name).

## Suggested Minimal First Milestone (1 PR)

- Add `src/geodata/mask_xarray.py` with:
  - `XarrayMask.from_name(...)`
  - `compute_grid_area(...)`
  - `attach(...)`
  - `apply(...)` (`where` + `multiply`)
- Reuse copied helper logic from `cutout.py` initially (refactor later).
- Add tests:
  - parity test with `Cutout.mask(...)` behavior for `attach()`
  - correctness test for `apply(...)`
  - area calculation sanity test

## Example Future Usage

```python
import geodata

# model output
ds_model = model.run(...)  # xr.Dataset with dims time/lat/lon (or x/y)

# load and align mask to ds_model grid
xmask = geodata.XarrayMask.from_name("china", grid=ds_model)

# 1) feature-style output (raw + mask + area)
attached = xmask.attach(ds_model, include_area=True)

# 2) direct masked output
masked = xmask.apply(ds_model, mode="where", include_area=True)
china_masked = masked["merged_mask"]
```

## Recommended Naming

- Keep existing `Mask` name for geospatial mask construction object.
- Use a distinct name for xarray adapter to avoid confusion:
  - preferred: `XarrayMask`
  - alternatives: `MaskApplier`, `MaskDatasetAdapter`

This separation keeps responsibilities clear:
- `Mask`: build/store masks
- `XarrayMask`: align/apply masks to model outputs
