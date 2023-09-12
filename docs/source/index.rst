.. Geodata documentation master file, created by
   sphinx-quickstart on Tue Aug 22 14:57:10 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to Geodata's documentation!
===================================

.. include:: intro.rst

.. toctree::
   :maxdepth: 1
   :caption: Quickstart
   :hidden:

   quick_start/packagesetup
   quick_start/input_output

.. toctree::
   :caption: Dataset Specific Tutorials
   :maxdepth: 1
   :hidden:

   datasets/era5/index
   datasets/merra2/index

.. toctree::
   :maxdepth: 1
   :caption: Mask
   :glob:
   :hidden:

   mask/*

.. toctree::
   :maxdepth: 1
   :caption: Parameterization
   :glob:
   :hidden:

   parameterizations/*

.. toctree::
   :maxdepth: 1
   :caption: Visualization
   :glob:
   :hidden:

   visualization/*

.. toctree::
   :maxdepth: 2
   :caption: API Reference
   :glob:
   :hidden:

   api_reference/*


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
