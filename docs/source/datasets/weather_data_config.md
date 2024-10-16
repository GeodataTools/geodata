# Advanced: Structure of Weather Data Configuration

In Geodata, every downloadable dataset are associated with a unique `(module, weather_data_config)` pair. 
In this tuple, the `module` typicallly refers to the source of dataset, while the `weather_data_config` is 
a dictionary that contains the information needed to download the specific form of the dataset. 

As Geodata currently supports `ERA5` and `MERRA2` modules, you can find all relevant weather data configuration
in each module's introduction pages here ([ERA5](era5/index.md), [MERRA2](merra2/index.md)). To find each config's actual definition, you can go to `src/geodata/datasets`. Within it, all available weather data configurations are located at the bottom of the file. 

In this tutorial, we will discuss the structure of each `weather_data_config` in more details.

## Example Weather Data Configuration

Below is an excerpt from MERRA2's `slv_flux_hourly` configuration:

```python
"slv_flux_hourly": dict(
    api_func=api_merra2,
    file_granularity="daily_multiple",
    tasks_func=tasks_daily_merra2,
    meta_prepare_func=prepare_meta_merra2,
    prepare_func=prepare_month_surface_flux,
    template=os.path.join(merra2_dir, "{year}/{month:0>2}/MERRA2_*.tavg1_2d_slv_flx_Nx.*.nc4"),
    url=[
        "https://goldsmr4.gesdisc.eosdis.nasa.gov/data/MERRA2/M2T1NXFLX.5.12.4/{year}/{month:0>2}/MERRA2_{spinup}.tavg1_2d_flx_Nx.{year}{month:0>2}{day:0>2}.nc4",
        "https://goldsmr4.gesdisc.eosdis.nasa.gov/data/MERRA2/M2T1NXSLV.5.12.4/{year}/{month:0>2}/MERRA2_{spinup}.tavg1_2d_slv_Nx.{year}{month:0>2}{day:0>2}.nc4",
    ],
    url_opendap=[
        "https://goldsmr4.gesdisc.eosdis.nasa.gov/opendap/MERRA2/M2T1NXFLX.5.12.4/{year}/{month:0>2}/MERRA2_{spinup}.tavg1_2d_flx_Nx.{year}{month:0>2}{day:0>2}.nc4",
        "https://goldsmr4.gesdisc.eosdis.nasa.gov/opendap/MERRA2/M2T1NXSLV.5.12.4/{year}/{month:0>2}/MERRA2_{spinup}.tavg1_2d_slv_Nx.{year}{month:0>2}{day:0>2}.nc4",
    ],
    fn=os.path.join(
        merra2_dir,
        "{year}/{month:0>2}/MERRA2_{spinup}.tavg1_2d_slv_flx_Nx.{year}{month:0>2}{day:0>2}.nc4",
    ),
    variables=[
        "ustar",
        "z0m",
        "disph",
        "rhoa",
        "ulml",
        "vlml",
        "tstar",
        "hlml",
        "tlml",
        "pblh",
        "hflux",
        "eflux",
        "u2m",
        "v2m",
        "u10m",
        "v10m",
        "u50m",
        "v50m",
    ],
    variables_list=[
        [
            "ustar",
            "z0m",
            "disph",
            "rhoa",
            "ulml",
            "vlml",
            "tstar",
            "hlml",
            "tlml",
            "pblh",
            "hflux",
            "eflux",
        ],
        ["u2m", "v2m", "u10m", "v10m", "u50m", "v50m"],
    ],
)
```

## Explanation of Each Field

### `api_func`

For each weather data configuration, we need to specify the API function that will be used to download the data. In our `slv_flux_hourly` example, we used the `api_merra2` in `merra2.py` to download the data. This function will be called by the `get_data` method of the dataset object with the following arguments `(toDownload, bounds, download_vars, product, product_type, downloadedFiles)`. If you encounter any issues with you try to download data from a specific module, this function will be a good place to debug. 

Typically, most weather data configurations within a module will share the same API function. However, there exists some exceptions, such as ERA5's `wind_3d_hourly`. Therefore, it's always a goodo idea to double-check which API function is being used for your desired weather data configuration.

### `file_granularity`

Different weather data configurations often may drastically differ in their "unit sizes". That is, for one month of data, some configurations may only take up 50 MiB of storage while others may take up 500 MiB. To avoid data corruption, Geodata stores data by periods. For example, if a given weather data configuration has a `daily_multiple` file granularity, then Geodata will store the data in daily "chunks" of xarray files. All supported file granularities are listed below:

- `daily` and `dailymeans`: Each file contains data for a single day.
- `monthly`: Each file contains data for a single month.
- `daily_multiple`: Each file contains data for a single day. However, each day's data maybe stored in multiple files. This is typical if a dataset is composed of data downloaded from multiple source datasets.
- `monthly_multiple`: Each file contains data for a single month. However, each month's data maybe stored in multiple files. This is typical if a dataset is composed of data downloaded from multiple source datasets.

### `url`

This field contains a list of URLs that will be used to download the data. The URLs are formatted using Python's [f-string](https://realpython.com/python-f-strings/) syntax.

### `url_opendap`

This field contains a list of URLs that will be used to download the data via [OPeNDAP](https://www.opendap.org/). The URLs are formatted using Python's [f-string](https://realpython.com/python-f-strings/) syntax.

### `fn`

This field represents the save path for individual files after downloading and processing. The path is formatted using Python's [f-string](https://realpython.com/python-f-strings/) syntax. It is worth noting that even if a weather data configuration has `multiple` file granularity (e.g. `daily_multiple`), the `fn` field will still only contain a single file path. This is because Geodata will automatically merge data from multiple sources into one after downloading and processing.

### `variables`

This field contains a list of variables in the xarray dataset that will be kept after download. 

#### Difference Between `variables` and `keywords`

Although not present here, some weather data configurations may also have a `keywords` field. This field aims to reconcile the difference between the variable names listed in the download request and the variable names in the downloaded dataset. For example, in ERA5 datasets, if a user requests the variable `100m_u_component_of_wind`, the downloaded dataset will contain a variable called `u100`. Thus, we need to specify `100m_u_component_of_wind` in the `keywords` field and `u100` in the `variables` field.

If the request variable names coincides with the name in the downloaded xarray files, then the `keywords` field is not needed. This is the case for MERRA2 datasets. Thus, you will not see a `keywords` field in MERRA2 weather data configurations.

### `variables_list`

This field is only used for weather data configurations with `multiple` file granularity (e.g. `daily_multiple`). It contains a list of lists of variables. Each inner list represents the variables that will be kept in a single file. For example, in the `slv_flux_hourly` example, the first inner list contains variables that will be kept in the first `slv` file, while the second inner list contains variables that will be kept in the second `flx` file.

<!-- Do we need to include documentation for all the unused fields? (e.g. `tasks_func`, `meta_prepare_func`, `prepare_func`, `template`) -->
