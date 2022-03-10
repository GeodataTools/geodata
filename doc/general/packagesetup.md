Back to the [Table of Contents](https://github.com/east-winds/geodata/blob/master/doc/general/tableofcontents.md).

# Geodata Package Setup

This guide covers how to install and configure the **geodata** package for local and cloud use.

Make sure that you have the following **required** software set up:

* [Python 3](https://www.python.org/downloads/)
* [Conda](https://docs.conda.io/projects/conda/en/latest/) (miniconda or anaconda) or [Pip](https://pip.pypa.io/en/stable/installation/) Python package management system.


Table of Contents:
- [Downloading and Configuring Geodata](#downloading-and-configuring-geodata)
- [Building Geodata](#building-geodata)
  - [Recommended Conda Installation Guide](#the-recommended-way-installation-in-conda-environment)
  - [MAC OS Pip Guide](#mac-os-installation-with-pip)
  - [Windows Pip Guide](#windows-installation-with-pip)
    - [pipwin instruction](#option-i-pipwin)
    - [wheel file instruction](#option-ii-wheel-files)

## Downloading and Configuring Geodata

To download **geodata**, open a terminal/shell window navigate to your preferred working directory, and run the following: (If you do not have Git installed, you may also directly download it with this [link](https://github.com/east-winds/geodata/archive/refs/heads/master.zip).

```
git clone https://github.com/east-winds/geodata.git
cd geodata
```

Before building the package, you'll first need to tell it two things: where to put/look for downloaded data, and where to store _cutouts_ - subsets of downloaded data needed to generate output variables, and _masks_ - geospatial raster files for suitability analysis.  To do so, copy `geodata/config-default.py` to a new file `geodata/config.py`.

* To configure where to store cutouts and masks, change the value of `cutout_dir` so that it points to a folder in your working directory.  You'll need to follow a format that fits your operating system:

```
## Mac OS
cutout_dir = '/Users/johndoe/desktop/geodata/data/cutouts'
mask_dir = '/Users/johndoe/desktop/geodata/data/masks'

## Windows
cutout_dir = 'C:/Users/johndoe/Desktop/geodata/data/cutouts'
mask_dir = 'C:/Users/johndoe/Desktop/geodata/data/masks'
```

**Note**: Make sure you are referencing a folder and path that already exist - the package currently does not create it for you.

* To configure where downloaded earth system data will be stored, change each datasets respective directory variable like so:

```
## Mac OS
merra2_dir = '/Users/johndoe/desktop/geodata/data/merra2'
era5_dir = '/Users/johndoe/desktop/geodata/data/era5'

## Windows
merra2_dir = 'C:/Users/johndoe/desktop/geodata/data/merra2'
era5_dir = 'C:/Users/johndoe/desktop/geodata/data/era5'
```

**Note**: Again, make sure you are referencing folders and paths that you've created beforehand - the package currently does not create them for you.

## Building Geodata

This guide includes 3 sections: the Conda environment (for Linux, Window, and Mac) guide, MAC OS with Pip guide, and Windows OS with Pip guide.

### The Recommended Way: Installation in Conda Environment

Conda is a powerful package manager and environment manager for Windows, macOS or Linux, and it provides easy installation for all operating systems. It is especially convenient if you are building Geodata on the cloud with potential installation permission issues.

If you already have Conda installed on your machine, jump straight to the `conda activate` step. Otherwise, you have 2 [options](https://conda.io/projects/conda/en/latest/user-guide/install/download.html#anaconda-or-miniconda): download Anaconda or miniconda. Installing Anaconda requires >3GB disk space and takes minutes to download, so we will choose **miniconda** instead because is a small, bootstrap version of Anaconda that includes only conda, Python, the packages they depend on, and a small number of other useful packages. You can get miniconda following use this [link](https://docs.conda.io/en/latest/miniconda.html#installing).

If you have conda 4.6 or later versions, in the terminal/shell, run the following command below to activate the conda [environment](https://docs.conda.io/projects/conda/en/latest/user-guide/getting-started.html#managing-environments).

```
conda activate
```

To use **geodata** in Python scripts by calling `import geodata`, you'll need to build the package.  To do so, in the terminal/shell window, navigate to the package's root directory (ie, "geodata"), and run the following:

```
python3 setup.py install
```

**Note I.**: If your ternimal prompts that python3 was not found, you can also try `python` or `py` instead and make sure to try `python --version` check for version.

**Note II.**: You will need to rebuild the package after making any changes to `config.py` in order for your changes to take effect.

**Note III.**: If running `python3 setup.py install` generates errors related to being unable to install the **rasterio** package due to conflicts with incompatible packages, you may need to reinstall Conda (either Anaconda or miniconda depending on what you went with during setup).  Then run the following commands:

```
conda update --all
```

```
conda install rasterio
```

### MAC OS Installation with Pip

In the MAC OS terminal, navigate to the package's root directory, and run the following:

```
python3 setup.py install
```

**Note I.**: If your ternimal prompts that python3 was not found, you can also try `python` or `py` instead and make sure to try `python --version` check for version.

**Note II.**: You will need to rebuild the package after making any changes to `config.py` in order for your changes to take effect.

The downside of not using a package management environment like conda is that we might run into dependency installation errors and have to fix it manuelly. (For a full list of dependencies, see the [setup.py](/setup.py) file or run the following command in the top directory: `python setup.py requirements`.) All dependencies should install automatically upon building the package, with possible exceptions such as the **rasterio** library, which requires Cython and GDAL. For the source of these instructions and more documentation about **rasterio**, see [the rasterio documentation](https://rasterio.readthedocs.io/en/latest/installation.html).

If one of the dependency, such as **rasterio** does not install automatically (we know this through the error message from the command above), we will have to install it seperately in the terminal:

```
pip install rasterio
```

Once the library above is successfully installed, re-run the installation command above to build Geodata:

```
python3 setup.py install
```

If there is an error message regarding one of Geodata's dependency, repeat the process and use `pip install` to seperately download it.

### Windows installation with Pip

In the Windows command prompt, navigate to the package's root directory, and run the following:

```
python3 setup.py install
```

**Note I.**: If your ternimal prompts that python3 was not found, you can also try `python` or `py` instead and make sure to try `python --version` check for version.

**Note II.**: You will need to rebuild the package after making any changes to `config.py` in order for your changes to take effect.

Similar to the MAC OS installation guide, it is likely that we would run into dependency installation error. We can figure out which dependency causes the error in the error message, and download that package seperately. For a full list of dependencies, see the [setup.py](/setup.py) file or run the following command in the top directory: `python setup.py requirements`.  All dependencies should install automatically upon building the package, with possible exceptions such as the **rasterio**  library, which requires other dependencies.

If one of the dependency, such as **GDAL** does not install automatically (we know this through the error message from the command above), we will have to install it seperately in the terminal. There are 2 options to solve this issue. Once we download the required dependency successfully, we can proceed by re-run the install command:

```
python3 setup.py install
```
#### Option I. Pipwin

This option is recommended. When we download libraries that is built on **GDAL**, we might run into this [issue](https://stackoverflow.com/q/54734667), where a GDAL API version must be specified.

**Pipwin** is a complementary tool for **pip** on Windows. **pipwin** installs unofficial python package binaries for windows provided by Christoph Gohlke [here](http://www.lfd.uci.edu/~gohlke/pythonlibs/).

Run the following commands to download **pipwin** and acquire the dependencies: (This solution is adopted from [stackoverflow](https://stackoverflow.com/a/58943939))

```
pip install pipwin
pipwin install shapely 
pipwin install gdal 
pipwin install fiona 
pipwin install pyproj 
pipwin install six 
pipwin install rtree 
```

You may need to also install the **wheel** package `pip install wheel` to facilitate building the wheels.

Similarly, if you run into installation errors regarding the **rasterio** or **bottleneck** packages, you can also call **pipwin install rasterio** or **pipwin install bottleneck** to download them.

#### Option II. Wheel Files

To install **rasterio** and the necessary GDAL library, we can download the appropriate binaries for your system by hand ([rasterio](https://www.lfd.uci.edu/~gohlke/pythonlibs/#rasterio) and [GDAL](https://www.lfd.uci.edu/~gohlke/pythonlibs/#gdal)) , place them into the current working directory, and run the following command in the downloads folder:

```
pip install -U pip
pip install {GDAL binary name here}.whl
pip install {rasterio binary name here}.whl
```

You may also need to also install the **wheel** package `pip install wheel` to facilitate building the wheels.

For more documentation/troubleshooting conda installation issues, see: https://github.com/conda-forge/rasterio-feedstock
