# `merge_layer` known issues (historical)

```{note}
**Historical context.** Older geodata versions could raise ``RasterioIOError: No such
file or directory`` when merging **in-memory** (``/vsimem``) mask layers. Current code
pins memory files for the lifetime of each layer reader so ``filter_layer`` →
``merge_layer`` normally works without pre-saving layers to disk.
```

## Symptom

``merge_layer`` (or ``merge_and`` / ``merge_sum`` after filters) fails with an error
referring to a missing path under ``/vsimem/``.

## Cause (legacy behavior)

Raster layers stored in GDAL memory files were sometimes closed before merge read them
back, so the virtual path was no longer valid.

## Current behavior

The mask module keeps layer readers alive while a ``Mask`` object uses in-memory
layers. If you still see this error on an old install, upgrade geodata or save
intermediate layers to disk before merging.

## Workaround (older versions)

1. Save filtered layers to GeoTIFF before ``merge_layer``.
2. Call ``save_mask(close_files=True)`` when finished, and avoid two ``Mask`` objects
   opening the same files simultaneously (see [mask troubleshooting](mask_troubleshoot.md)).

## Tests

Regression coverage lives under ``tests/pr/mask/test_mask_merge_inmemory.py``.
