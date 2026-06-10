Legacy workflow
===============

The pages below document the original Geodata API built around
``Dataset``, ``Cutout``, ``geodata.convert``, and Cutout-based masking,
including MERRA2 download and cutout tutorials.

.. note::

   This path is **not** part of the current tested workflow
   (``load_dataset`` → models → ``XarrayMask``). It remains available for
   existing analyses and reference.

For the recommended path, see the :doc:`documentation homepage </index>`.

.. toctree::
   :maxdepth: 1

   workflow
   mask_on_cutout
   merra2/index
   merra2/merra2_download
   merra2/merra2_outputs
   merra2/merra2
   wind_extrapolation
