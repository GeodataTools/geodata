# Copyright 2020 Michael Davidson (UCSD), William Honaker.

# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License as
# published by the Free Software Foundation; either version 3 of the
# License, or (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.


"""
GEODATA

Geospatial Data Collection and "Pre-Analysis" Tools
"""

import glob
import os
from calendar import monthrange
from tempfile import mkstemp

import numpy as np
import requests
import xarray as xr
from requests.exceptions import HTTPError
from tqdm.contrib.logging import tqdm_logging_redirect

from ..config import merra2_dir
from ..logging import logger

datadir = merra2_dir

# Model and Projection Settings
projection = "latlong"


def convert_and_subset_lons_lats_merra2(ds, xs, ys):
    # Rename geographic dimensions to x,y
    # Subset x,y according to xs, ys

    if not isinstance(xs, slice):
        first, second, last = np.asarray(xs)[[0, 1, -1]]
        xs = slice(first - 0.1 * (second - first), last + 0.1 * (second - first))
    if not isinstance(ys, slice):
        first, second, last = np.asarray(ys)[[0, 1, -1]]
        ys = slice(first - 0.1 * (second - first), last + 0.1 * (second - first))

    ds = ds.sel(lat=ys)

    # Lons should go from -180. to +180.
    if len(ds.coords["lon"].sel(lon=slice(xs.start + 360.0, xs.stop + 360.0))):
        ds = xr.concat(
            [ds.sel(lon=slice(xs.start + 360.0, xs.stop + 360.0)), ds.sel(lon=xs)],
            dim="lon",
        )
        ds = ds.assign_coords(
            lon=np.where(
                ds.coords["lon"].values <= 180,
                ds.coords["lon"].values,
                ds.coords["lon"].values - 360.0,
            )
        )
    else:
        ds = ds.sel(lon=xs)

    ds = ds.rename({"lon": "x", "lat": "y"})
    ds = ds.assign_coords(lon=ds.coords["x"], lat=ds.coords["y"])
    return ds


def subset_x_y_merra2(ds, xs, ys):
    # Subset x,y according to xs, ys
    # Assumes convert_and_subset_lons_lats_merra2 already run

    if not isinstance(xs, slice):
        first, second, last = np.asarray(xs)[[0, 1, -1]]
        xs = slice(first - 0.1 * (second - first), last + 0.1 * (second - first))
    if not isinstance(ys, slice):
        first, second, last = np.asarray(ys)[[0, 1, -1]]
        ys = slice(first - 0.1 * (second - first), last + 0.1 * (second - first))

    ds = ds.sel(y=ys)
    ds = ds.sel(x=xs)

    return ds


def _rename_and_clean_coords(ds, add_lon_lat=True):
    """Rename 'longitude' and 'latitude' columns to 'x' and 'y'

    Optionally (add_lon_lat, default:True) preserves latitude and longitude columns as 'lat' and 'lon'.
    """

    ds = ds.rename({"lon": "x", "lat": "y"})
    if add_lon_lat:
        ds = ds.assign_coords(lon=ds.coords["x"], lat=ds.coords["y"])
    return ds


def api_merra2(toDownload, fileGranularity, downloadedFiles):
    if len(toDownload) == 0:
        logger.info("All MERRA2 files for this dataset have been downloaded.")
    else:
        multi = bool(fileGranularity in ("daily_multiple", "monthly_multiple"))

        total = len(toDownload) * (len(toDownload[0]) - 2) if multi else len(toDownload)
        with tqdm_logging_redirect(
            loggers=[logger], total=total, dynamic_ncols=True
        ) as pbar:
            for f in toDownload:
                error_files = []
                if multi:
                    fd, target = mkstemp(suffix=".nc4")
                else:
                    fd = 0
                    target = f[1]
                os.makedirs(os.path.dirname(f[1]), exist_ok=True)
                logger.info("Preparing API calls for %s", f[1])
                logger.info("Making request to %s", f[2])
                result = requests.get(f[2], timeout=30)
                try:
                    result.raise_for_status()
                    with open(target, "wb") as fout:
                        fout.write(result.content)
                except HTTPError as http_err:
                    logger.warning("HTTP error occurred: %s", http_err)  # Python 3.6
                    error_files.append(f[1])
                except Exception as err:  # pylint: disable=broad-except
                    logger.warning("Other error occurred: %s", err)
                    error_files.append(f[1])
                pbar.update()

                if multi:
                    # ds_main = xr.open_dataset(target)

                    temp_files = [[fd, target]]
                    for k in range(3, len(f)):
                        logger.info("Making request to %s", f[k])
                        result = requests.get(f[k], timeout=30)
                        fd_temp, target_temp = mkstemp(suffix=".nc4")
                        temp_files.append([fd_temp, target_temp])
                        try:
                            result.raise_for_status()
                            with open(target_temp, "wb") as fout:
                                fout.write(result.content)
                        except HTTPError as http_err:
                            logger.warning(
                                "HTTP error occurred: %s", http_err
                            )  # Python 3.6
                            error_files.append(f[k])
                        except Exception as err:  # pylint: disable=broad-except
                            logger.warning("Other error occurred: %s", err)
                            error_files.append(f[k])
                        # ds_toadd = xr.open_dataset(target_temp)
                        # ds_main = xr.merge([ds_main, ds_toadd])
                        # os.close(fd_temp)
                        # os.unlink(target_temp)
                        pbar.update()

                    ds_main = xr.open_mfdataset(
                        [fn[1] for fn in temp_files], combine="by_coords"
                    )
                    ds_main.to_netcdf(f[1])
                    ds_main.close()
                    # ds_toadd.close()  # close last xr open file

                    # close and clear temp files
                    # os.close(fd)
                    # os.unlink(target)
                    for tf in temp_files:
                        os.close(tf[0])
                        os.unlink(tf[1])

                if len(error_files) > 0:
                    logger.warning("Unsuccessful download for %s", error_files)
                else:
                    logger.info("Successfully downloaded data for %s", f[1])
                    downloadedFiles.append((f[0], f[1]))


def prepare_meta_merra2(xs, ys, year, month, template, module, **params):
    # 	Load dataset into metadata

    # fn = next(glob.iglob(template.format(year=year, month=month)))
    # with xr.open_dataset(fn) as ds:
    # 	ds = ds.coords.to_dataset()
    # 	ds = convert_and_subset_lons_lats_merra2(ds, xs, ys)
    # 	meta = ds.load()

    # Set spinup variable (see MERRA2 documentation, p. 13)
    spinup = spinup_year(year, month)

    fns = glob.iglob(template.format(year=year, month=month, spinup=spinup))
    with xr.open_mfdataset(fns, combine="by_coords") as ds:
        ds = ds.coords.to_dataset()
        ds = convert_and_subset_lons_lats_merra2(ds, xs, ys)
        meta = ds.load()

    return meta


def prepare_month_surface_flux(fn, year, month, xs, ys):
    if not os.path.isfile(fn):
        return None
    with xr.open_dataset(fn) as ds:
        logger.info("Opening %s", fn)
        # logger.info("Cutout dims: %s", ds.dims)
        # logger.info("Cutout coords: %s", ds.coords)

        ds = _rename_and_clean_coords(ds)

        # logger.info("Cutout dims: %s", ds.dims)
        # logger.info("Cutout coords: %s", ds.coords)

        ds = subset_x_y_merra2(ds, xs, ys)

        # logger.info("Cutout dims: %s", ds.dims)
        # logger.info("Cutout coords: %s", ds.coords)

        # some variable renaming
        try:
            # 	z0m=roughness
            # 	wind variables not in wndXXm format
            ds = ds.rename({"z0m": "roughness"})
        except Exception as e:
            logger.warning("Unable to rename variables in %s. Exception: %s", fn, e)

        ds["wndlml"] = np.sqrt(ds["ulml"] ** 2 + ds["vlml"] ** 2).assign_attrs(
            units=ds["ulml"].attrs["units"], long_name="LML wind speed"
        )
        if "tlml" in list(ds.data_vars):
            ds["temperature"] = ds["tlml"]

        yield (year, month), ds


def prepare_month_aerosol(fn, year, month, xs, ys):
    if not os.path.isfile(fn):
        return None
    with xr.open_dataset(fn) as ds:
        logger.info("Opening %s", fn)
        ds = _rename_and_clean_coords(ds)
        ds = subset_x_y_merra2(ds, xs, ys)
        yield (year, month), ds


def prepare_dailymeans_surface_flux(fn, year, month, xs, ys):
    if not os.path.isfile(fn):
        return None
    with xr.open_dataset(fn) as ds:
        logger.info("Opening %s", fn)
        # logger.info("Cutout dims: %s", ds.dims)
        # logger.info("Cutout coords: %s", ds.coords)

        ds = _rename_and_clean_coords(ds)

        # logger.info("Cutout dims: %s", ds.dims)
        # logger.info("Cutout coords: %s", ds.coords)

        ds = subset_x_y_merra2(ds, xs, ys)

        # logger.info("Cutout dims: %s", ds.dims)
        # logger.info("Cutout coords: %s", ds.coords)

        # some variable renaming
        try:
            # 	z0m=roughness
            # 	wind variables not in wndXXm format
            ds = ds.rename({"t2mmean": "temperature", "tprecmax": "precipitation"})
        except Exception as e:
            logger.warning("Unable to rename variables in %s. Exception: %s", fn, e)

        # ['HOURNORAIN', 'T2MMAX', 'T2MMEAN', 'T2MMIN', 'TPRECMAX']

        yield (year, month), ds


def prepare_slv_radiation(fn, year, month, xs, ys):
    if not os.path.isfile(fn):
        return None
    with xr.open_dataset(fn) as ds:
        logger.info("Opening %s", fn)
        # logger.info("Cutout dims: %s", ds.dims)
        # logger.info("Cutout coords: %s", ds.coords)

        ds = _rename_and_clean_coords(ds)

        # logger.info("Cutout dims: %s", ds.dims)
        # logger.info("Cutout coords: %s", ds.coords)

        ds = subset_x_y_merra2(ds, xs, ys)

        # logger.info("Cutout dims: %s", ds.dims)
        # logger.info("Cutout coords: %s", ds.coords)
        try:
            ds = ds.rename(
                {
                    "albedo": "albedo",
                    "swgdn": "influx",
                    "swtdn": "influx_toa",
                    "t2m": "temperature",
                }
            )
        except Exception as e:
            logger.warning("Unable to rename variables in %s. Exception: %s", fn, e)
        yield (year, month), ds


## TODO def prepare_month_radiation
# with np.errstate(divide='ignore', invalid='ignore'):
# 	ds['albedo'] = (((ds['ssrd'] - ds['ssr'])/ds['ssrd']).fillna(0.)
# 					.assign_attrs(units='(0 - 1)', long_name='Albedo'))
# ds['influx_diffuse'] = ((ds['ssrd'] - ds['influx_direct'])
# 						.assign_attrs(units='J m**-2',
# 									long_name='Surface diffuse solar radiation downwards'))
# ds = ds.drop(['ssrd', 'ssr'])
#
# # Convert from energy to power J m**-2 -> W m**-2 and clip negative fluxes
# for a in ('influx_direct', 'influx_diffuse', 'influx_toa'):
# 	ds[a] = ds[a].clip(min=0.) / (60.*60.)
# 	ds[a].attrs['units'] = 'W m**-2'


def tasks_daily_merra2(xs, ys, yearmonths, prepare_func, **meta_attrs):
    if not isinstance(xs, slice):
        xs = slice(*xs.values[[0, -1]])
    if not isinstance(ys, slice):
        ys = slice(*ys.values[[0, -1]])
    fn = meta_attrs["fn"]

    logger.info(yearmonths)
    logger.info(
        [
            (year, month, day)
            for year, month in yearmonths
            for day in range(1, monthrange(year, month)[1] + 1, 1)
        ]
    )

    return [
        dict(
            prepare_func=prepare_func,
            xs=xs,
            ys=ys,
            year=year,
            month=month,
            fn=fn.format(
                year=year, month=month, day=day, spinup=spinup_year(year, month)
            ),
        )
        for year, month in yearmonths
        for day in range(1, monthrange(year, month)[1] + 1, 1)
    ]


def tasks_monthly_merra2(xs, ys, yearmonths, prepare_func, **meta_attrs):
    if not isinstance(xs, slice):
        xs = slice(*xs.values[[0, -1]])
    if not isinstance(ys, slice):
        ys = slice(*ys.values[[0, -1]])
    fn = meta_attrs["fn"]

    logger.info(yearmonths)
    logger.info([(year, month) for year, month in yearmonths])

    return [
        dict(
            prepare_func=prepare_func,
            xs=xs,
            ys=ys,
            year=year,
            month=month,
            fn=fn.format(year=year, month=month, spinup=spinup_year(year, month)),
        )
        for year, month in yearmonths
    ]


weather_data_config = {
    "surface_flux_hourly": dict(
        api_func=api_merra2,
        file_granularity="daily",
        tasks_func=tasks_daily_merra2,
        meta_prepare_func=prepare_meta_merra2,
        prepare_func=prepare_month_surface_flux,
        template=os.path.join(
            merra2_dir, "{year}/{month:0>2}/MERRA2_*.tavg1_2d_flx_Nx.*.nc4"
        ),
        url="https://goldsmr4.gesdisc.eosdis.nasa.gov/data/MERRA2/M2T1NXFLX.5.12.4/{year}/{month:0>2}/MERRA2_{spinup}.tavg1_2d_flx_Nx.{year}{month:0>2}{day:0>2}.nc4",
        url_opendap="https://goldsmr4.gesdisc.eosdis.nasa.gov/opendap/MERRA2/M2T1NXFLX.5.12.4/{year}/{month:0>2}/MERRA2_{spinup}.tavg1_2d_flx_Nx.{year}{month:0>2}{day:0>2}.nc4.nc4",
        fn=os.path.join(
            merra2_dir,
            "{year}/{month:0>2}/MERRA2_{spinup}.tavg1_2d_flx_Nx.{year}{month:0>2}{day:0>2}.nc4",
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
        ],
    ),
    "slv_flux_hourly": dict(
        api_func=api_merra2,
        file_granularity="daily_multiple",
        tasks_func=tasks_daily_merra2,
        meta_prepare_func=prepare_meta_merra2,
        prepare_func=prepare_month_surface_flux,
        template=os.path.join(
            merra2_dir, "{year}/{month:0>2}/MERRA2_*.tavg1_2d_slv_flx_Nx.*.nc4"
        ),
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
    ),
    "surface_flux_monthly": dict(
        api_func=api_merra2,
        file_granularity="monthly",
        tasks_func=tasks_monthly_merra2,
        meta_prepare_func=prepare_meta_merra2,
        prepare_func=prepare_month_surface_flux,
        template=os.path.join(merra2_dir, "{year}/MERRA2_*.tavgM_2d_flx_Nx.*.nc4"),
        url="https://goldsmr4.gesdisc.eosdis.nasa.gov/data/MERRA2_MONTHLY/M2TMNXFLX.5.12.4/{year}/MERRA2_{spinup}.tavgM_2d_flx_Nx.{year}{month:0>2}.nc4",
        fn=os.path.join(
            merra2_dir, "{year}/MERRA2_{spinup}.tavgM_2d_flx_Nx.{year}{month:0>2}.nc4"
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
        ],
    ),
    "surface_flux_dailymeans": dict(
        api_func=api_merra2,
        file_granularity="dailymeans",
        tasks_func=tasks_daily_merra2,
        meta_prepare_func=prepare_meta_merra2,
        prepare_func=prepare_dailymeans_surface_flux,
        template=os.path.join(
            merra2_dir, "{year}/{month:0>2}/MERRA2_*.statD_2d_slv_Nx.*.nc4"
        ),
        url="https://goldsmr4.gesdisc.eosdis.nasa.gov/data/MERRA2/M2SDNXSLV.5.12.4/{year}/{month:0>2}/MERRA2_{spinup}.statD_2d_slv_Nx.{year}{month:0>2}{day:0>2}.nc4",
        fn=os.path.join(
            merra2_dir,
            "{year}/{month:0>2}/MERRA2_{spinup}.statD_2d_slv_Nx.{year}{month:0>2}{day:0>2}.nc4",
        ),
        variables=["hournorain", "tprecmax", "t2mmax", "t2mmean", "t2mmin"],
    ),
    "slv_radiation_hourly": dict(
        api_func=api_merra2,
        file_granularity="daily_multiple",
        tasks_func=tasks_daily_merra2,
        meta_prepare_func=prepare_meta_merra2,
        prepare_func=prepare_slv_radiation,
        template=os.path.join(
            merra2_dir, "{year}/{month:0>2}/MERRA2_*.tavg1_2d_slv_rad_Nx.*.nc4"
        ),
        url=[
            "https://goldsmr4.gesdisc.eosdis.nasa.gov/data/MERRA2/M2T1NXSLV.5.12.4/{year}/{month:0>2}/MERRA2_{spinup}.tavg1_2d_slv_Nx.{year}{month:0>2}{day:0>2}.nc4",
            "https://goldsmr4.gesdisc.eosdis.nasa.gov/data/MERRA2/M2T1NXRAD.5.12.4/{year}/{month:0>2}/MERRA2_{spinup}.tavg1_2d_rad_Nx.{year}{month:0>2}{day:0>2}.nc4",
        ],
        url_opendap=[
            "https://goldsmr4.gesdisc.eosdis.nasa.gov/opendap/MERRA2/M2T1NXSLV.5.12.4/{year}/{month:0>2}/MERRA2_{spinup}.tavg1_2d_slv_Nx.{year}{month:0>2}{day:0>2}.nc4.nc4",
            "https://goldsmr4.gesdisc.eosdis.nasa.gov/opendap/MERRA2/M2T1NXRAD.5.12.4/{year}/{month:0>2}/MERRA2_{spinup}.tavg1_2d_rad_Nx.{year}{month:0>2}{day:0>2}.nc4.nc4",
        ],
        fn=os.path.join(
            merra2_dir,
            "{year}/{month:0>2}/MERRA2_{spinup}.tavg1_2d_slv_rad_Nx.{year}{month:0>2}{day:0>2}.nc4",
        ),
        variables=["albedo", "swgdn", "swtdn", "t2m"],
        variables_list=[["t2m"], ["albedo", "swgdn", "swtdn"]],
    ),
    "slv_radiation_monthly": dict(
        api_func=api_merra2,
        file_granularity="monthly_multiple",
        tasks_func=tasks_monthly_merra2,
        meta_prepare_func=prepare_meta_merra2,
        prepare_func=prepare_slv_radiation,
        template=os.path.join(merra2_dir, "{year}/MERRA2_*.tavgM_2d_slv_rad_Nx.*.nc4"),
        url=[
            "https://goldsmr4.gesdisc.eosdis.nasa.gov/data/MERRA2_MONTHLY/M2TMNXSLV.5.12.4/{year}/MERRA2_{spinup}.tavgM_2d_slv_Nx.{year}{month:0>2}.nc4",
            "https://goldsmr4.gesdisc.eosdis.nasa.gov/data/MERRA2_MONTHLY/M2TMNXRAD.5.12.4/{year}/MERRA2_{spinup}.tavgM_2d_rad_Nx.{year}{month:0>2}.nc4",
        ],
        fn=os.path.join(
            merra2_dir,
            "{year}/MERRA2_{spinup}.tavgM_2d_slv_rad_Nx.{year}{month:0>2}.nc4",
        ),
        variables=["albedo", "swgdn", "swtdn", "t2m"],
    ),
    "surface_aerosol_hourly": dict(
        api_func=api_merra2,
        file_granularity="daily",
        tasks_func=tasks_daily_merra2,
        meta_prepare_func=prepare_meta_merra2,
        prepare_func=prepare_month_aerosol,
        template=os.path.join(
            merra2_dir, "{year}/{month:0>2}/MERRA2_*.tavg1_2d_aer_Nx.*.nc4"
        ),
        url="https://goldsmr4.gesdisc.eosdis.nasa.gov/data/MERRA2/M2T1NXAER.5.12.4/{year}/{month:0>2}/MERRA2_{spinup}.tavg1_2d_aer_Nx.{year}{month:0>2}{day:0>2}.nc4",
        fn=os.path.join(
            merra2_dir,
            "{year}/{month:0>2}/MERRA2_{spinup}.tavg1_2d_aer_Nx.{year}{month:0>2}{day:0>2}.nc4",
        ),
        variables=["bcsmass", "dusmass25", "ocsmass", "so4smass", "sssmass25"],
    ),
}

# os.path.join(merra2_dir, '{year}/MERRA2_{spinup}.tavgM_2d_flx_Nx.{year}{month:0>2}.nc4'),

# list of routines in weather_data_config to download wind data
wind_files = ["surface_flux"]

# TODO: same for solar
solar_files = []

## Whatever is calling this needs to be directed to the correct weather config instead
# meta_data_config = dict(prepare_func=prepare_meta_merra2,
# 						 template=os.path.join(merra2_dir, '{year}/{month:0>2}/MERRA2_*.tavg1_2d_flx_Nx.*.nc4'))

# Separate files for each day (coded in weather_data_config list)
# daily_files = True # needs to be specified somewhere else

# Latitude stored south to north (ie forward, = True) or north to south
lat_direction = True

# Spinup variable
spinup_var = True


def spinup_year(year, month):
    if year >= 1980 and year < 1992:
        spinup = "100"
    elif year >= 1992 and year < 2001:
        spinup = "200"
    elif year >= 2001 and year < 2011:
        spinup = "300"
    elif year >= 2011 and year < 2020:
        spinup = "400"
    elif year == 2020 and month == 9:
        spinup = "401"
    else:
        spinup = "400"

    return spinup
