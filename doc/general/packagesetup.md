# Geodata Package Setup

How to install and configure the **geodata** package for local use.
This guide assumes you have the following installed/setup:

* [Python 3](https://www.python.org/downloads/)
* [Git](https://git-scm.com/downloads) installed and a [Github](https://github.com/) account.
* (Recommended) A development environment for writing and running Python scripts and/or Jupyter notebooks. A couple good options are: [Visual Studio Code](https://code.visualstudio.com/) with [Jupyter Notebook integration](https://code.visualstudio.com/docs/python/jupyter-support), or [JupyterLab](https://jupyterlab.readthedocs.io/en/stable/).

## Installing and Configuring Geodata

To install **geodata**, open a terminal/shell window navigate to your preferred working directory, and run the following:

```
git clone https://github.com/east-winds/geodata.git
cd geodata
```

Before building the package, you'll first need to tell it three things: where to put/look for downloaded data, where to store _cutouts_--subsets of downloaded data needed to generate output variables--and where to stores _masks_--geospatial layers and manipulations to extract on cutouts.  To do so, copy `geodata/config-default.py` to a new file `geodata/config.py`, and adjust the following lines.

* To configure where to store cutouts, change the value of `cutout_dir` so that it points to a folder in your working directory.  You'll need to follow a format that fits your operating system:

Mac OS

```
cutout_dir = '/Users/johndoe/desktop/geodata/data/cutouts'
```

Windows

```
cutout_dir = 'C:/Users/johndoe/Desktop/geodata/data/cutouts'
```

(**Note**: Make sure you are referencing a folder and path that already exists.)

* To configure where to store masks, change the value of `mask_dir` so that it points to a folder in your working directory.

* To configure where downloaded earth system data will be stored, change each datasets respective directory variable like so:

For MERRA2:

```
## Mac OS
merra2_dir = '/Users/johndoe/desktop/geodata/data/merra2'

## Windows
merra2_dir = 'C:/Users/johndoe/desktop/geodata/data/merra2'

```


## Building Geodata

To use **geodata**, you'll need to build the package.  To do so, open a terminal/shell window, navigate to the package's root directory (ie, "geodata"), and run the following:

```
python3 setup.py install
```

This will build the package and allow you to use it in Python scripts by calling `import geodata`.

**Note**: You will need to rebuild the package after making any changes to `config.py` in order for your changes to take effect.

## Dependencies and GDAL Driver Installation

For a full list of dependencies, see the [setup.py](/setup.py) file or run the following command in the top directory: `python setup.py requirements`.  All dependencies should install automatically upon building the package, with the possible exception of **rasterio**, which requires GDAL.  If **rasterio** does not install automatically, follow the instructions below specific to your operating system + Python package manager.  For the source of these instructions and more documentation about **rasterio**, see [the rasterio documentation](https://rasterio.readthedocs.io/en/latest/installation.html).

### Mac OS with pip

To install **rasterio** and the necessary GDAL library, simply run:

```
pip install rasterio
```

### Windows with pip

To install **rasterio** and the necessary GDAL library, first download the appropriate binaries for your system ([rasterio](https://www.lfd.uci.edu/~gohlke/pythonlibs/#rasterio) and [GDAL](https://www.lfd.uci.edu/~gohlke/pythonlibs/#gdal)) and run the following command in the downloads folder:

```
pip install -U pip
pip install {GDAL binary name here}.whl
pip install {rasterio binary name here}.whl
```

You may need to also install the **wheel** package `pip install wheel` to facilitate building the wheels.

### Anaconda

To install, enable the `conda-forge` channel and install.

```
conda config --add channels conda-forge
conda install rasterio
```

For more documentation/troubleshooting conda installation issues, see: https://github.com/conda-forge/rasterio-feedstock
