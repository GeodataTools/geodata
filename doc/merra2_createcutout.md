# Creating Cutouts with MERRA2 data

A short guide on how to use **geodata** to create cutout - subsets of data based on specific time and geographic ranges - for MERRA2 data.

## Setup

To start, import the required dependencies:

```import logging
logging.basicConfig(level=logging.INFO)
import atlite
```

`import atlite` is required to use **geodata**, while launching a logger allows for detailed debugging via the console.

## Preparing the cutout

A cutout is the basis for any data or analysis output by the **geodata** package.  Cutouts are stored in the directory `cutout_dir` configured in `config.py` (to set up `config.py`, see here LINK HERE)

```cutout = atlite.Cutout(name="merra2-europe-sub24-2011-01",
                       module="merra2",
                       xs=slice(30, 41.56244222),
                       ys=slice(33.56459975, 35),
                       years=slice(2011, 2011),
                       months=slice(1,1))
```

To prepare a cutout, the following must be specified for `atlite.Cutout()`:

* The cutout name
* The source dataset
* Time range
* Geographic range as represented by `xy` coordinates.

The example in the code block above uses MERRA2 data, as specified by the `module` parameter.

```module="merra2"
```

`xs` and `ys` in combination with the `slice()` function allow us to specify a geographic range based on longitude and latitude.  The above example subsets a portion of Europe.

`years` and `months` are used to subset the time range.  For both functions, the first value represents the start point, and the second value represents the end point.  The above example creates a cutout for January 2011.


## Creating a Cutout

To create a cutout, run:
```cutout.prepare();
```
The **geodata** package will create a folder in the cutout directory (`cutout_dir`) you specified in `config.py` with the name specified in `atlite.Cutout()` (in the above example, `merra2-europe-sub24-2011-01`).  The folder, depending on the date range, will then contain one or more netcdf files containing data from the original MERRA2 files falling within the temporal and geographical ranges indicated when the cutout was created.

Before running `cutout.prepare()`, **geodata** will check for the presence of the cutout and abort if the cutout already exists.  If you want to force the regeneration of a cutout, run the command with the parameter `overwrite=True`.


TODO: Add info about examining cutout data here.