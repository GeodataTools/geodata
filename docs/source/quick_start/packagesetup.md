# Installation

This guide covers how to install and configure the **geodata** package for local and cloud use.

Make sure that you have the following software set up:

* [Python 3](https://www.python.org/downloads/)

* Package Management System (Optional but highly recommended to isolate individual projects dependent packages)
  - [conda](https://docs.conda.io/projects/conda/en/latest/) (miniconda or anaconda)


```{note}

The majority of the content below applies to macOS and Linux. Any possible differences for Windows users are noted.
Nontheless, it is worth noting that using **geodata** on Windows may present some additional challenges due to its
file access properties.

Thus, it is recommended to use geodata a Linux-based system. If you are using Windows, you may
want to consider using the [Windows Subsystem for Linux (WSL)](https://learn.microsoft.com/en-us/windows/wsl/install)
or a virtual machine.
```

## Downloading Geodata

To install **geodata**, open a terminal/shell window and run the following:

```bash
pip install geodata-re
```

### Advanced: Installing from Source

If you want to install **geodata** from source, you can clone the repository and install it using the following commands:

```bash
git clone https://github.com/GeodataTools/geodata.git
pip install geodata/
```

## Configuring File Storage Location

To configure where to store downloaded and processed files, define an environment variable called `GEODATA_ROOT` and save in your shell configuration files, such as `.bashrc` or `.zshrc`:

```bash
export GEODATA_ROOT=<YOUR_PATH_HERE>
```

```{note}
If you are using a Windows machine, you can set the environment variable by running the following command in the command prompt:

```bash
setx GEODATA_ROOT <YOUR PATH HERE>
```

If you are running geodata in a Jupyter Notebook, you could define the variable by adding and running the following cell:
```
%setenv GEODATA_ROOT <YOUR PATH HERE>
```

If you do not define this variable, all datasets and cutouts will be stored under `~/.local/geodata` by default.

[ADD typical folders to be created under GEODATA_ROOT]

## Anaconda/miniconda Environment

[Anaconda](https://www.anaconda.com/download)/[miniconda](https://docs.conda.io/en/latest/miniconda.html) is a powerful package manager and environment manager for Windows, macOS or Linux, and it provides easy installation for all operating systems. It is especially convenient if you are building Geodata on the cloud with potential installation permission issues.

If you already have Anaconda/miniconda installed on your machine, jump straight to the `conda activate` step. Otherwise, you have 2 [options](https://conda.io/projects/conda/en/latest/user-guide/install/download.html#anaconda-or-miniconda): download Anaconda or miniconda. Installing Anaconda requires >3GB disk space and takes minutes to download, so we will choose **miniconda** instead because is a small, bootstrap version of Anaconda that includes only conda, Python, the packages they depend on, and a small number of other useful packages.

If you have conda 4.6 or later versions, in the terminal/shell, run the following command below to activate the conda [environment](https://docs.conda.io/projects/conda/en/latest/user-guide/getting-started.html#managing-environments).

```bash

conda create --name geodata python=3.12
conda activate <ENVIRONMENT_NAME>
pip install geodata-re
```

Once you activate the environment, any packages you install (including geodata) will be isolated.
If you have a new project and wish to start over again, you can create a new environment and install the package again.

```bash
conda deactivate
conda env remove --name geodata
```
