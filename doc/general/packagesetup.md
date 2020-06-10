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

Before building the package, you'll first need to tell it two things: where to put/look for downloaded data, and where to store _cutouts_--subsets of downloaded data needed to generate output variables.  To do so, copy `geodata/config-default.py` to a new file `geodata/config.py`.

* To configure where to store cutouts, change the value of `cutout_dir` so that it points to a folder in your working directory like so:
```
cutout_dir = '/Users/johndoe/desktop/geodata/data/cutouts'
```
(**Note**: Make sure you are referencing a folder and path that already exists - the package currently does not create it for you.)

* To configure where downloaded earth system data will be stored, change each datasets respective directory variable like so:

For MERRA2:
```
merra2_dir = '/Users/johndoe/desktop/geodata/data/merra2'
```
For ERA5:
```
era5_dir = '/Users/johndoe/desktop/geodata/data/era5'
```
(**Note**: Again, make sure you are referencing folders and paths that you've created beforehand - the package currently does not create them for you.)

## Building Geodata
To use **geodata**, you'll need to build the package.  To do so, open a terminal/shell window, navigate to the package's root directory (ie, "geodata"), and run the following:

```
python3 setup.py install
```

This will build the package and allow you to use it in Python scripts by calling `import geodata`.

For a full list of requirements, see the [setup.py](/setup.py) file or run the following command in the top directory: `python setup.py requirements`.

**Note**: You will need to rebuild the package after making any changes to `config.py` in order for your changes to take effect.