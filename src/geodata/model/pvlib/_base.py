# Copyright 2016-2017 Gorm Andresen (Aarhus University), Jonas Hoersch (FIAS), Tom Brown (FIAS)
# Copyright 2020 Michael Davidson (UCSD), William Honaker, Jiahe Feng (UCSD), Yuanbo Shi
# Copyright 2023-2024 Xiqiang Liu, 2025 Keyu Long

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

TODO: Documentation here

"""
import os
import platform
import pandas as pd
import xarray as xr
import time
from multiprocessing import Pool, Manager, cpu_count as mp_cpu_count
from pvlib import pvsystem
from pvlib.location import Location
from pvlib.modelchain import ModelChain
from timezonefinder import TimezoneFinder

from .._base import BaseModel, _normalize_slice_for_sel, _should_use_parallel_reading
from geodata.logging import logger
from .calculations import calculate_pvlib_solarposition, calculate_ghi, calculate_relative_humidity, calculate_precipitable_water, convert_kelvin_to_celsius
from tqdm.auto import tqdm


class ModelChainConfig:
    """
    Defines pvlib ModelChain parameters as a class that
    can be passed to one or more instances of pvlib_model().
    Allows user to reuse a common set of ModelChain parameters across multiple
    PVSystems or even multiple cutouts.  

    Parameters
    ----------
    clearsky_model : string, default 'ineichen'
        Specifies the clear-sky model. Passed to location.get_clearsky. 
        Only used when DNI is not found in the weather inputs.
    transposition_model : string, default 'haydavies'
        Specifies the transposition model. Passed to system.get_irradiance.
    solar_position_method : string, default 'nrel_numpy'
        Specifies the method for calculating solar positions. Passed to location.get_solarposition.
    airmass_model : string, default 'kastenyoung1989'
        Specifies the airmass model. Passed to location.get_airmass.
    dc_model : string or function, optional
        Specifies the DC model. Valid strings are 'sapm', 'desoto', 'cec', 'pvsyst', 'pvwatts'. 
        If not specified, the model will be inferred from the parameters of system.arrays[i].module_parameters. 
        A user-defined function may also be provided, with the ModelChain instance passed as the first argument.
    ac_model : string or function, optional
        Specifies the AC model. Valid strings are 'sandia', 'adr', 'pvwatts'. 
        If not specified, the model will be inferred from the parameters of system.inverter_parameters. 
        A user-defined function may also be provided, with the ModelChain instance passed as the first argument.
    aoi_model : string or function, optional
        Specifies the angle of incidence (AOI) model. Valid strings are 'physical', 'ashrae', 'sapm', 'martin_ruiz', 
        'interp', 'no_loss'. If not specified, the model will be inferred from the parameters of 
        system.arrays[i].module_parameters. A user-defined function may also be provided, 
        with the ModelChain instance passed as the first argument.
    spectral_model : string or function, optional
        Specifies the spectral model. Valid strings are 'sapm', 'first_solar', 'no_loss'. 
        If not specified, the model will be inferred from the parameters of system.arrays[i].module_parameters. 
        A user-defined function may also be provided, with the ModelChain instance passed as the first argument.
    temperature_model : string or function, optional
        Specifies the temperature model. Valid strings are 'sapm', 'pvsyst', 'faiman', 'fuentes', 'noct_sam'. 
        A user-defined function may also be provided, with the ModelChain instance passed as the first argument.
    dc_ohmic_model : string or function, default 'no_loss'
        Specifies the DC ohmic loss model. Valid strings are 'dc_ohms_from_percent', 'no_loss'. 
        A user-defined function may also be provided, with the ModelChain instance passed as the first argument.
    losses_model : string or function, default 'no_loss'
        Specifies the losses model. Valid strings are 'pvwatts', 'no_loss'. 
        A user-defined function may also be provided, with the ModelChain instance passed as the first argument.
    name : string, optional
        Specifies the name of the ModelChain instance.
    
    For full documentation, see:
    - pvlib.modelchain.ModelChain(): 
        https://pvlib-python.readthedocs.io/en/stable/reference/generated/pvlib.modelchain.ModelChain.html

    """
    def __init__(
        self,
        clearsky_model='ineichen',
        transposition_model='haydavies',
        solar_position_method='nrel_numpy',
        airmass_model='kastenyoung1989',
        dc_model=None,
        ac_model=None,
        aoi_model=None,
        spectral_model=None,
        temperature_model=None,
        dc_ohmic_model='no_loss',
        losses_model='no_loss',
        name=None
    ):
        self.clearsky_model = clearsky_model
        self.transposition_model = transposition_model
        self.solar_position_method = solar_position_method
        self.airmass_model = airmass_model
        self.dc_model = dc_model
        self.ac_model = ac_model
        self.aoi_model = aoi_model
        self.spectral_model = spectral_model
        self.temperature_model = temperature_model
        self.dc_ohmic_model = dc_ohmic_model
        self.losses_model = losses_model
        self.name = name

    def model_chain_to_kwargs(self):
        return self.__dict__


def _detect_available_cpus() -> int:
    """
    Detect the number of available CPUs using multiple methods.
    
    Tries multiple detection methods in order:
    1. SLURM environment variables (if in SLURM job) - authoritative for HPC clusters
    2. psutil (if available) - most reliable, respects CPU affinity
    3. Linux cgroups v2 (if available) - respects container limits
    4. Linux cgroups v1 (if available) - respects container limits
    5. multiprocessing.cpu_count() - standard library fallback
    6. os.cpu_count() - last resort
    7. Defaults to 1 if all methods fail
    
    Returns:
        int: Number of available CPUs (at least 1)
    """
    detected_cpus = None
    method_used = None
    
    # Method 1: Try SLURM environment variables (for HPC clusters)
    # SLURM is authoritative when present, so check this first
    try:
        # Check if we're in a SLURM job
        if os.getenv("SLURM_JOB_ID") is not None:
            # Try SLURM_CPUS_PER_TASK first (most common and reliable)
            slurm_cpus_per_task = os.getenv("SLURM_CPUS_PER_TASK")
            if slurm_cpus_per_task is not None:
                try:
                    detected_cpus = int(slurm_cpus_per_task)
                    if detected_cpus > 0:
                        method_used = "SLURM_CPUS_PER_TASK"
                        logger.debug(
                            f"_detect_available_cpus: Detected {detected_cpus} CPU(s) "
                            f"using SLURM_CPUS_PER_TASK={slurm_cpus_per_task}"
                        )
                        # SLURM is authoritative - return immediately
                        logger.debug(
                            f"_detect_available_cpus: Successfully detected {detected_cpus} CPU(s) "
                            f"using method: {method_used}"
                        )
                        return detected_cpus
                except (ValueError, TypeError):
                    logger.debug(
                        f"_detect_available_cpus: SLURM_CPUS_PER_TASK={slurm_cpus_per_task} "
                        f"is not a valid integer, trying other SLURM variables"
                    )
            
            # If SLURM_CPUS_PER_TASK not available, try SLURM_JOB_CPUS_PER_NODE
            if detected_cpus is None:
                slurm_job_cpus = os.getenv("SLURM_JOB_CPUS_PER_NODE")
                if slurm_job_cpus is not None:
                    try:
                        # SLURM_JOB_CPUS_PER_NODE can be a comma-separated list for multi-node jobs
                        # Take the first value (current node)
                        cpus_str = slurm_job_cpus.split(',')[0]
                        detected_cpus = int(cpus_str)
                        if detected_cpus > 0:
                            method_used = "SLURM_JOB_CPUS_PER_NODE"
                            logger.debug(
                                f"_detect_available_cpus: Detected {detected_cpus} CPU(s) "
                                f"using SLURM_JOB_CPUS_PER_NODE={slurm_job_cpus} "
                                f"(using first node value)"
                            )
                            # SLURM is authoritative - return immediately
                            logger.debug(
                                f"_detect_available_cpus: Successfully detected {detected_cpus} CPU(s) "
                                f"using method: {method_used}"
                            )
                            return detected_cpus
                    except (ValueError, TypeError, IndexError):
                        logger.debug(
                            f"_detect_available_cpus: SLURM_JOB_CPUS_PER_NODE={slurm_job_cpus} "
                            f"could not be parsed, trying other methods"
                        )
            
            # If still not found, try SLURM_CPUS_ON_NODE (but this is less reliable)
            # as it shows CPUs on node, not necessarily allocated to job
            if detected_cpus is None:
                slurm_cpus_on_node = os.getenv("SLURM_CPUS_ON_NODE")
                if slurm_cpus_on_node is not None:
                    try:
                        # Can be a comma-separated list for multi-node jobs
                        cpus_str = slurm_cpus_on_node.split(',')[0]
                        detected_cpus = int(cpus_str)
                        if detected_cpus > 0:
                            method_used = "SLURM_CPUS_ON_NODE"
                            logger.debug(
                                f"_detect_available_cpus: Detected {detected_cpus} CPU(s) "
                                f"using SLURM_CPUS_ON_NODE={slurm_cpus_on_node} "
                                f"(using first node value). Note: This may not reflect "
                                f"actual CPU allocation to the job."
                            )
                            # SLURM is authoritative - return immediately
                            logger.debug(
                                f"_detect_available_cpus: Successfully detected {detected_cpus} CPU(s) "
                                f"using method: {method_used}"
                            )
                            return detected_cpus
                    except (ValueError, TypeError, IndexError):
                        logger.debug(
                            f"_detect_available_cpus: SLURM_CPUS_ON_NODE={slurm_cpus_on_node} "
                            f"could not be parsed, trying other methods"
                        )
            
            if detected_cpus is None:
                logger.debug(
                    "_detect_available_cpus: Running in SLURM job (SLURM_JOB_ID present) "
                    "but no usable CPU count variables found. Trying other detection methods."
                )
    except Exception as e:
        logger.debug(f"_detect_available_cpus: SLURM detection failed: {e}, trying other methods")
    
    # Method 2: Try psutil (most reliable, respects CPU affinity and cgroups)
    if detected_cpus is None:
        try:
            import psutil
            detected_cpus = psutil.cpu_count(logical=False)  # Physical cores first
            if detected_cpus is None or detected_cpus == 0:
                detected_cpus = psutil.cpu_count(logical=True)  # Fallback to logical cores
            if detected_cpus is not None and detected_cpus > 0:
                method_used = "psutil"
                logger.debug(f"_detect_available_cpus: Detected {detected_cpus} CPU(s) using psutil")
        except ImportError:
            logger.debug("_detect_available_cpus: psutil not available, trying other methods")
        except Exception as e:
            logger.debug(f"_detect_available_cpus: psutil failed: {e}, trying other methods")
    
    # Method 3: Try Linux cgroups v2 (for containers)
    if detected_cpus is None and platform.system() == "Linux":
        try:
            # Check cgroup v2 cpu.max (format: "max" or "quota period")
            cgroup_path = "/sys/fs/cgroup/cpu.max"
            if os.path.exists(cgroup_path):
                with open(cgroup_path, 'r') as f:
                    content = f.read().strip()
                    if content != "max":
                        parts = content.split()
                        if len(parts) == 2:
                            quota = int(parts[0])
                            period = int(parts[1])
                            if quota > 0 and period > 0:
                                detected_cpus = max(1, int(quota / period))
                                method_used = "cgroups_v2"
                                logger.debug(
                                    f"_detect_available_cpus: Detected {detected_cpus} CPU(s) "
                                    f"using cgroups v2 (quota={quota}, period={period})"
                                )
        except Exception as e:
            logger.debug(f"_detect_available_cpus: cgroups v2 check failed: {e}")
        
        # Method 4: Try Linux cgroups v1 (for containers)
        if detected_cpus is None:
            try:
                # Check cgroup v1 cpu.cfs_quota_us and cpu.cfs_period_us
                quota_path = "/sys/fs/cgroup/cpu/cpu.cfs_quota_us"
                period_path = "/sys/fs/cgroup/cpu/cpu.cfs_period_us"
                if os.path.exists(quota_path) and os.path.exists(period_path):
                    with open(quota_path, 'r') as f:
                        quota = int(f.read().strip())
                    with open(period_path, 'r') as f:
                        period = int(f.read().strip())
                    if quota > 0 and period > 0:
                        detected_cpus = max(1, int(quota / period))
                        method_used = "cgroups_v1"
                        logger.debug(
                            f"_detect_available_cpus: Detected {detected_cpus} CPU(s) "
                            f"using cgroups v1 (quota={quota}, period={period})"
                        )
            except Exception as e:
                logger.debug(f"_detect_available_cpus: cgroups v1 check failed: {e}")
    
    # Method 5: Try multiprocessing.cpu_count()
    if detected_cpus is None:
        try:
            detected_cpus = mp_cpu_count()
            if detected_cpus is not None and detected_cpus > 0:
                method_used = "multiprocessing.cpu_count()"
                logger.debug(
                    f"_detect_available_cpus: Detected {detected_cpus} CPU(s) "
                    f"using multiprocessing.cpu_count()"
                )
        except Exception as e:
            logger.debug(f"_detect_available_cpus: multiprocessing.cpu_count() failed: {e}")
    
    # Method 6: Try os.cpu_count() as last resort
    if detected_cpus is None:
        try:
            detected_cpus = os.cpu_count()
            if detected_cpus is not None and detected_cpus > 0:
                method_used = "os.cpu_count()"
                logger.debug(
                    f"_detect_available_cpus: Detected {detected_cpus} CPU(s) using os.cpu_count()"
                )
        except Exception as e:
            logger.debug(f"_detect_available_cpus: os.cpu_count() failed: {e}")
    
    # Final fallback: default to 1
    if detected_cpus is None or detected_cpus <= 0:
        detected_cpus = 1
        method_used = "default_fallback"
        logger.debug(
            "_detect_available_cpus: All detection methods failed. "
            "Defaulting to 1 CPU and logging debug message."
        )
        logger.debug(
            "_detect_available_cpus: This may indicate the program is running in a restricted "
            "environment (container, cgroup limits, or CPU affinity restrictions)."
        )
    else:
        logger.debug(
            f"_detect_available_cpus: Successfully detected {detected_cpus} CPU(s) "
            f"using method: {method_used}"
        )
    
    return detected_cpus


def _process_single_coordinate(args):
    """
    Helper function to process a single coordinate for multiprocessing.
    
    This function must be at module level to be picklable for multiprocessing.
    
    Args:
        args: Tuple containing:
            - coord: Tuple of (y, x) coordinates
            - weather_data: DataFrame with weather data
            - system: PVSystem object
            - model_chain_kwargs: Dictionary of ModelChain configuration
            - ptc: Module PTC value
            - n_mods: Number of modules per string
            - progress_dict: Shared dictionary for progress tracking (optional)
            - coord_index: Index of this coordinate in the total list
            - total_coords: Total number of coordinates to process
            - compact_output: Whether to keep only ac/pv in output
    
    Returns:
        Tuple of (coord, subset_df) where subset_df contains the processed data
    """
    (
        (y, x),
        weather_data,
        system,
        model_chain_kwargs,
        ptc,
        n_mods,
        progress_dict,
        coord_index,
        total_coords,
        compact_output,
    ) = args
    
    try:
        # Extract subset for this coordinate
        subset = weather_data.loc[(slice(None), y, x), :].reset_index(['x', 'y'])
        
        # Get timezone
        tz_str = TimezoneFinder().timezone_at(lat=y, lng=x)
        if tz_str is None:
            raise ValueError(f"Timezone not found for coordinates ({y}, {x})")
        
        # Create location
        location = Location(latitude=y, longitude=x, tz=tz_str)  # type: ignore[arg-type]
        
        # Create and run ModelChain
        mc = ModelChain(
            system,
            location,
            **model_chain_kwargs
        )
        mc.run_model(subset)
        
        # Calculate outputs
        subset['ac'] = mc.results.ac
        subset.loc[subset['ac'] < 0, 'ac'] = 0
        subset['pv'] = subset['ac'] / (ptc * n_mods)
        
        # Update progress if progress_dict is provided
        if progress_dict is not None:
            with progress_dict['lock']:
                progress_dict['completed'] += 1
                completed = progress_dict['completed']
                elapsed = time.time() - progress_dict['start_time']
                
                # Log progress periodically
                log_interval = max(1, min(100, total_coords // 10))
                if completed % log_interval == 0 or completed == total_coords:
                    avg_time_per_coord = elapsed / completed if completed > 0 else 0
                    remaining_coords = total_coords - completed
                    eta = avg_time_per_coord * remaining_coords
                    progress_dict['last_log'] = {
                        'completed': completed,
                        'total': total_coords,
                        'coord': (y, x),
                        'elapsed': elapsed,
                        'avg_time': avg_time_per_coord,
                        'eta': eta
                    }
                    progress_dict['should_log'] = True

        # Re-pack results into a MultiIndex so that
        # xr.Dataset.from_dataframe() reconstructs x and y as dimensions.
        #
        # `subset` is currently indexed only by `time` (x/y were reset into columns),
        # which would otherwise cause the output to have only `time` as a coordinate.
        if compact_output:
            subset_out = subset[['ac', 'pv']].copy()
        else:
            subset_out = subset.copy()
        subset_out = subset_out.assign(y=y, x=x)
        subset_out = subset_out.reset_index()

        # After reset_index(), the time column name can vary (e.g. 'time' vs 'index').
        if subset.index.name is None:
            subset_out = subset_out.rename(columns={'index': 'time'})
        elif subset.index.name != 'time':
            subset_out = subset_out.rename(columns={subset.index.name: 'time'})

        subset_out = subset_out.set_index(['time', 'x', 'y'])
        return (y, x), subset_out
        
    except Exception as e:
        logger.error(f"Error processing coordinate ({y}, {x}): {str(e)}")
        raise


class Pvlib(BaseModel):
    """The pvlib model"""
    @property
    def type(self) -> str:
        return "pvlib"

    SUPPORTED_WEATHER_DATA_CONFIGS = ("wind_solar_hourly", "wind_solar_hourly_test")

    @property
    def prepared(self) -> bool:
        """This model does not need to be prepared"""
        return True
    
    def prepare(self, force: bool = False):
        """Skip preparation - this model doesn't need it."""
        logger.info("This model does not require preparation. Skipping.")
        return
    
    def init_model_config(
            self,
            clearsky_model='ineichen',
            transposition_model='haydavies',
            solar_position_method='nrel_numpy',
            airmass_model='kastenyoung1989',
            dc_model=None,
            ac_model=None,
            aoi_model=None,
            spectral_model=None,
            temperature_model=None,
            dc_ohmic_model='no_loss',
            losses_model='no_loss',
            name=None
        ):
            self.config = ModelChainConfig(
                clearsky_model= clearsky_model,
                transposition_model= transposition_model, 
                solar_position_method= solar_position_method,
                airmass_model= airmass_model,
                dc_model= dc_model,
                ac_model= ac_model, 
                aoi_model= aoi_model,
                spectral_model= spectral_model,
                temperature_model= temperature_model,
                dc_ohmic_model= dc_ohmic_model,
                losses_model= losses_model,
                name= name
            )

    def retrieve_sam(self, samfile, path=None):
        """
        Wrapper for pvlib.pvsystem.retrieve_sam(). Retrieves latest module 
        and inverter info from a file bundled with pvlib, a path or a 
        URL (like SAM’s website), and returns it as a Pandas DataFrame.

        Supported databases:
        - CEC module database
        - Sandia Module database
        - CEC Inverter database
        - Anton Driesse Inverter database

        Parameters
        ----------
        name : string
                Use one of the following strings to retrieve a database bundled with pvlib:
                    - ’CECMod’ - returns the CEC module database
                    - ’CECInverter’ - returns the CEC Inverter database
                    - ’SandiaInverter’ - returns the CEC Inverter database 
                        (CEC is only current inverter db available; tag kept for backwards compatibility)
                    - ’SandiaMod’ - returns the Sandia Module database
                    - ’ADRInverter’ - returns the ADR Inverter database
        
        Optional Parameters
        ----------
        path : string
                Path to a CSV file or a URL.

        Returns: DataFrame

        See also:
            - pvlib.pvsystem.retrieve_sam(): 
                https://pvlib-python.readthedocs.io/en/stable/reference/generated/pvlib.pvsystem.retrieve_sam.html

        """
        return pvsystem.retrieve_sam(name=samfile, path=path)
    
    def init_pv_system(self, *args, **kwargs):
        """
        Wrapper for pvlib.pvsystem.PVSystem().
        The PVSystem class defines a standard set of PV system attributes
        and modeling functions. This class describes the collection and 
        interactions of PV system components rather than an installed system
        on the ground. It is typically used in combination with Location 
        and ModelChain objects.

        The class supports basic system topologies consisting of:
            - N total modules arranged in series (modules_per_string=N, strings_per_inverter=1).
            - M total modules arranged in parallel (modules_per_string=1, strings_per_inverter=M).
            - NxM total modules arranged in M strings of N modules each 
            (modules_per_string=N, strings_per_inverter=M).

        For full documentation, see: https://pvlib-python.readthedocs.io/en/stable/reference/generated/pvlib.pvsystem.PVSystem.html

        Parameters
        ----------
        arrays : array (optional)
                An Array or list of arrays that are part of the system. 
                See pvlib documentation for full description.
        surface_tilt : float
                Surface tilt angles in decimal degrees. The tilt angle is 
                defined as degrees from horizontal (e.g. surface facing up = 0, 
                surface facing horizon = 90).
        surface_azimuth : float
                Azimuth angle of the module surface. North=0, East=90, South=180, West=270.
        albedo : float
                Ground surface albedo. If not supplied, then surface_type is used to look up 
                a value in pvlib.albedo.SURFACE_ALBEDOS. If surface_type is also not supplied 
                then a ground surface albedo of 0.25 is used.
        surface_type : string
                The ground surface type. See pvlib.albedo.SURFACE_ALBEDOS for valid values.
        module : string
                The model name of the modules. May be used to look up the module_parameters dictionary via some other method.
        module_type : string 
                Describes the module’s construction. Valid strings are ‘glass_polymer’ and ‘glass_glass’. 
                Used for cell and module temperature calculations.
        module_parameters : dict
                Module parameters as defined by the SAPM, CEC, or other.
        temperature_model_parameters : dict
                Temperature model parameters as required by one of the models in pvlib.temperature (excluding poa_global, temp_air and wind_speed).
        modules_per_string : int, float
                See system topology discussion above.
        strings_per_inverter : int, float 
                See system topology discussion above.
        inverter : string 
                The model name of the inverters. May be used to look up the inverter_parameters dictionary via some other method.
        inverter_parameters : dict
                Inverter parameters as defined by the SAPM, CEC, or other.
        racking_model : string 
                Valid strings are ‘open_rack’, ‘close_mount’, and ‘insulated_back’. 
                Used to identify a parameter set for the SAPM cell temperature model.
        losses_parameters : dict 
                Losses parameters as defined by PVWatts or other.    
        name : string (optional)
        
        """
        self.pv_system = pvsystem.PVSystem(*args, **kwargs)


    def _estimate_dataset(self, params: xr.Dataset, **kwargs) -> xr.Dataset | xr.DataArray:  # type: ignore[override]
        """Estimate PV output from prepared dataset.
        
        Args:
            params: Dataset (already filtered by years/months/xs/ys from BaseModel)
            **kwargs: Additional parameters (not used currently, but available)
        
        Returns:
            Dataset with AC power and PV capacity (returns Dataset, but BaseModel expects DataArray)
        """

        compact_output = kwargs.get("compact_output", True)
        result = self._pvlib_model(
            params,
            self.pv_system,
            self.config,
            compact_output=compact_output,
        )
        return result
    
    def estimate(self,
        years: slice | None = None,
        months: slice | None = None,
        xs: slice | None = None,
        ys: slice | None = None,
        compact_output: bool = True,
        **kwargs,
        ) -> xr.DataArray:
            """Get pvlib model results.
            
            This method processes data month-by-month to avoid memory issues with large datasets.
            Results from each month are concatenated along the time dimension.
            
            Args:
                years: Year range (slice)
                months: Month range (slice)
                xs: X-coordinate range (slice)
                ys: Y-coordinate range (slice)
                compact_output: If True (default), return only `ac` and `pv`
                    as data variables. If False, keep full per-coordinate output.
                **kwargs: Additional parameters
            
            Returns:
                Dataset with AC power and PV capacity, concatenated across all months
            """
            if getattr(self, 'pv_system', None) is None:
                raise ValueError("pv_system is not initialized. Call init_pv_system() first.")
            if getattr(self, 'config', None) is None:
                raise ValueError("model_config is not initialized. Call init_model_config() first.")
            
            # Get result objects for the requested time range
            if years is None and months is None:
                results = self.flattened_results
            elif months is None:
                # If years specified but months not, use all months
                if years is None:
                    results = self.flattened_results
                else:
                    results = self.get_result_year_month(years, slice(1, 13))
            else:
                # Both years and months specified
                if years is None:
                    # If only months specified, need to get all years
                    # Use the source dataset's year range
                    years = self.source.years
                results = self.get_result_year_month(years, months)
            
            if not results:
                raise ValueError("No results found for the specified year/month range.")
            
            # Process month-by-month to manage memory
            logger.info(
                f"Processing {len(results)} month(s) month-by-month to manage memory usage"
            )
            
            engine = "h5netcdf"
            parallel = _should_use_parallel_reading()
            
            monthly_results = []
            
            for result in tqdm(results, desc="Processing months", unit="month"):
                # Load only this month's raw data files
                ref_files = result.ref_files
                
                if not ref_files:
                    logger.warning(
                        f"No files found for {result.year:04d}-{result.month:02d}, skipping."
                    )
                    continue
                
                logger.debug(
                    f"Loading {len(ref_files)} file(s) for {result.year:04d}-{result.month:02d} "
                    f"with engine={engine}, parallel={parallel}"
                )
                
                # Open this month's dataset
                with xr.open_mfdataset(
                    ref_files,
                    engine=engine,
                    parallel=parallel,
                ) as params:
                    # Rename coordinates to x/y if they use longitude/latitude naming
                    # This must happen before spatial filtering
                    # Check dimensions first (for .sel() to work), then coordinates
                    rename_dict = {}
                    if "longitude" in params.dims:
                        rename_dict["longitude"] = "x"
                    elif "lon" in params.dims:
                        rename_dict["lon"] = "x"
                    if "latitude" in params.dims:
                        rename_dict["latitude"] = "y"
                    elif "lat" in params.dims:
                        rename_dict["lat"] = "y"
                    
                    if rename_dict:
                        params = params.rename(rename_dict)
                    
                    # Apply spatial filtering if specified.
                    # Normalize slice order for descending coordinates (e.g. ERA5 latitude);
                    # otherwise .sel() returns empty.
                    if xs is not None:
                        x_slice = _normalize_slice_for_sel(params.coords["x"], xs) if "x" in params.coords else xs
                        params = params.sel(x=x_slice)
                    if ys is not None:
                        y_slice = _normalize_slice_for_sel(params.coords["y"], ys) if "y" in params.coords else ys
                        params = params.sel(y=y_slice)
                    
                    # Transform raw dataset to standardized format
                    # This applies the same transformations as prepare_func
                    # (renames variables, calculates derived quantities, etc.)
                    dataset_cls = type(self.source)
                    if hasattr(dataset_cls, 'transform_wind_solar_dataset'):
                        # Call the classmethod to transform the dataset
                        params = dataset_cls.transform_wind_solar_dataset(params)  # type: ignore[attr-defined]
                    else:
                        logger.warning(
                            "Dataset does not have transform_wind_solar_dataset method. "
                            "Assuming data is already in the correct format."
                        )
                    
                    # Process this month's data
                    monthly_output = self._estimate_dataset(
                        params,
                        compact_output=compact_output,
                        **kwargs,
                    )
                    
                    # Store the result (will concatenate later)
                    monthly_results.append(monthly_output)
            
            if not monthly_results:
                raise ValueError("No data was successfully processed for the specified range.")
            
            # Concatenate all monthly results along the time dimension
            logger.info(f"Concatenating {len(monthly_results)} month(s) of results")
            
            # Ensure all datasets have compatible coordinates
            # Sort by time to ensure proper ordering
            combined_result = xr.concat(monthly_results, dim='time')
            
            # Sort by time to ensure chronological order
            if 'time' in combined_result.coords:
                combined_result = combined_result.sortby('time')
            
            # Standardize output dimension order across models:
            # `("time", "x", "y")`.
            desired_order = ("time", "x", "y")
            if all(d in combined_result.dims for d in desired_order):
                combined_result = combined_result.transpose(*desired_order)
            return combined_result
        
    def _prepare_pvlib_ds(self, ds: xr.Dataset, *varnames: str) -> xr.Dataset:
        """
        Prepares an `xarray.Dataset` from a geodata `cutout` class for use in model simulations using `pvlib`.
        This function extracts specified variables from the `cutout` dataset, calculates additional parameters 
        like global horizontal irradiance (GHI), precipitable water, and solar position, and renames fields to 
        align with expected inputs.
        
        PARALLELIZATION OPTIONS:
        ------------------------
        This function processes the entire dataset at once using vectorized operations (xarray/numpy).
        Most operations are already parallelized at the numpy level (BLAS/MKL).
        
        Potential optimizations:
        1. If dataset is very large, consider chunking by coordinates and processing in parallel
        2. The solar position calculation (calculate_pvlib_solarposition) could be parallelized
           across time steps if it's not already vectorized
        3. Consider using dask arrays for lazy evaluation if memory is a concern

        Requires a cutout with the following variables:

        - **influx_diffuse** (*float*) - Diffuse horizontal irradiance.  
        - **influx_direct** (*float*) - Direct normal irradiance.  
        - **dewpoint_temperature** (*float*) - Dewpoint temperature in Celsius.  
        - **temperature** (*float*) - Air temperature in Celsius.  
        - **wnd100m** (*float*) - Wind speed at 100m.  

        Outputs an `xarray.Dataset` with the following variables:

        - **dhi** (*float*) - Diffuse horizontal irradiance.  
        - **dni** (*float*) - Direct normal irradiance.  
        - **ghi** (*float*) - Global horizontal irradiance (calculated via :code:`_calculate_ghi()`).  
        - **temp_air** (*float*) - Air temperature in Celsius.  
        - **wind_speed** (*float*) - Wind speed at 100m.  
        - **precipitable_water** (*float*) - Precipitable water (calculated via :code:`_calculate_precipitable_water()`).

        Parameters
        ----------
        ds : xarray.Dataset
            Must contain following variables: influx_diffuse, influx_direct,
            dewpoint_temperature, temperature, wnd100m.
        varnames : string
            String values representing names of required variables.

        Returns
        -------
        weather_data : `xarray.Dataset`
            Dataset containing necessary variables to run `pvlib` model simulations.

        """
        prepare_start = time.time()
        
        if varnames:
            # Check which variables are actually available
            available_vars = [v for v in varnames if v in ds.data_vars]
            missing_vars = [v for v in varnames if v not in ds.data_vars]
            if missing_vars:
                logger.warning(f"Missing variables: {missing_vars}. Available: {list(ds.data_vars.keys())}")
            if available_vars:
                ds = ds[available_vars]
            else:
                logger.error(f"None of the requested variables {varnames} are available in dataset")
                raise KeyError(f"None of the requested variables {varnames} are available. Available variables: {list(ds.data_vars.keys())}")

        temperature_celsius = convert_kelvin_to_celsius(ds.temperature)

        relative_humidity = calculate_relative_humidity(
            temperature_celsius,
            convert_kelvin_to_celsius(ds.dewpoint_temperature),
            # convert_kelvin_to_celsius(ds.d2m),
        )

        precipitable_water = calculate_precipitable_water(
            temperature_celsius,
            relative_humidity
        )

        sp = calculate_pvlib_solarposition(ds)
        ghi = calculate_ghi(ds, sp['zenith'])

        ds = (
            ds
            .assign(
                ghi=ghi,
                temperature=temperature_celsius,
                precipitable_water=precipitable_water
            )
            .rename({
                'influx_diffuse': 'dhi',
                'influx_direct': 'dni',
                'temperature': 'temp_air',
                'wnd100m': 'wind_speed'
            })
        )

        result = ds[[
            "dhi", 
            "dni",
            "ghi", 
            "temp_air", 
            "wind_speed", 
            "precipitable_water"
        ]]
        
        prepare_time = time.time() - prepare_start
        logger.info(f"_prepare_pvlib_ds: {prepare_time:.2f}s")
        
        return result
    
    def _pvlib_model(
        self,
        ds: xr.Dataset, 
        system: pvsystem.PVSystem,
        model_chain_config: ModelChainConfig,
        vars: list[str] = ["influx_diffuse", "influx_direct", "dewpoint_temperature", "temperature", "wnd100m"],
        n_jobs: int | None = None,
        compact_output: bool = True,
    ) -> xr.Dataset:
        
        """
        Applies a `pvlib` model using :code:`pvlib.modelchain.ModelChain()` across all unique coordinates 
        represented in a `geodata` cutout. This function prepares input weather data, initializes the 
        `pvlib` model, and runs simulations for each set of coordinates, outputting an xarray dataset 
        containing all simulation results.
        
        PARALLELIZATION OPTIONS:
        ------------------------
        This function processes coordinates sequentially, which is the main bottleneck.
        Recommended parallelization approaches:
        
        Option 1: multiprocessing.Pool (Recommended for CPU-bound tasks)
        ----------------------------------------
        - Use multiprocessing.Pool to process coordinates in parallel
        - Each worker processes a subset of coordinates independently
        - Pros: True parallelism, good for CPU-bound pvlib calculations
        - Cons: Requires pickling system/model_chain_config objects, higher memory usage
        
        Option 2: concurrent.futures.ThreadPoolExecutor
        ------------------------------------------------
        - Use threads for I/O-bound operations (if any)
        - Less overhead than multiprocessing
        - Pros: Lower memory overhead, faster startup
        - Cons: Limited by GIL for CPU-bound tasks (pvlib is CPU-bound, so not ideal)
        
        Option 3: concurrent.futures.ProcessPoolExecutor
        ------------------------------------------------
        - Similar to multiprocessing.Pool but with a simpler API
        - Pros: Cleaner API, better error handling
        - Cons: Similar to Option 1
        
        Option 4: joblib.Parallel
        --------------------------
        - High-level parallel processing library
        - Pros: Simple API, good progress reporting, handles pickling well
        - Cons: Additional dependency
        
        Option 5: Dask (for distributed computing)
        -------------------------------------------
        - For very large datasets across multiple machines
        - Pros: Scales to clusters, handles memory efficiently
        - Cons: More complex setup, overhead for small datasets
        
        Implementation suggestion:
        - Create a helper function: _process_single_coordinate(y, x, weather_data, system, model_chain_config, ptc, n_mods)
        - Use multiprocessing.Pool.map() or ProcessPoolExecutor.map() to parallelize
        - Consider chunking coordinates into batches to balance load
        - Use n_jobs parameter to control parallelism (default: os.cpu_count())

        Requires a cutout with the following variables:

        - **influx_diffuse** (*float*) - Diffuse horizontal irradiance.  
        - **influx_direct** (*float*) - Direct normal irradiance.  
        - **dewpoint_temperature** (*float*) - Dewpoint temperature in Celsius.  
        - **temperature** (*float*) - Air temperature in Celsius.  
        - **wnd100m** (*float*) - Wind speed at 100m.  

        Outputs an `xarray.Dataset` containing:

        - **ac** (*float*) - AC photovoltaic output (W).
        - **pv** (*float*) - Photovoltaic capacity.

        Parameters
        ----------
        cutout : geodata **cutout** class
            Cutout generated by the `geodata` library, based on the ERA5 dataset.  
            Must contain the required meteorological variables.
        system : pvlib **PVSystem** class
            The photovoltaic system to be simulated.  Generated by :code:`geodata.pvlib.pv_system()`
        model_chain_config : `ModelChainConfig`
            Configuration object for :code:`pvlib.modelchain.ModelChain()` with model parameters.
        vars : list of str, optional
            List of variable names required for simulation. Defaults to:
            ['influx_diffuse', 'influx_direct', 'dewpoint_temperature', 'temperature', 'wnd100m'].
        compact_output : bool, optional
            If True (default), keep only `ac` and `pv` data variables in the
            final output. If False, keep the full per-coordinate output.

        Returns
        -------
        xr.Dataset
            Dataset containing ac power output and pv capacity across all coordinates in the cutout.

        """
        ptc = system.arrays[0].module_parameters['PTC']
        n_mods = system.arrays[0].modules_per_string

        weather_data = self._prepare_pvlib_ds(ds, *vars).to_dataframe()
        unique_coords = weather_data.index.droplevel('time').drop_duplicates()
        total_coords = len(unique_coords)
        
        # Determine number of workers using robust CPU detection
        if n_jobs is None:
            n_jobs = _detect_available_cpus()
            logger.debug(
                f"_pvlib_model: Auto-detected {n_jobs} available CPU(s) for parallel processing"
            )
        else:
            logger.debug(
                f"_pvlib_model: Using user-specified n_jobs={n_jobs} for parallel processing"
            )
        
        # Ensure n_jobs is valid: at least 1, and not more than total coordinates
        n_jobs = max(1, min(n_jobs, total_coords))
        
        if n_jobs == 1:
            logger.debug(
                "_pvlib_model: Using sequential processing (n_jobs=1). "
                "This may be due to: only 1 coordinate, CPU detection returned 1, "
                "or user specified n_jobs=1"
            )
        
        logger.info(
            f"Processing {total_coords} coordinate(s) using {n_jobs} worker process(es)"
        )
        
        # Prepare arguments for parallel processing
        model_chain_kwargs = model_chain_config.model_chain_to_kwargs()
        
        # Create shared progress tracking dictionary
        manager = Manager()
        progress_dict = manager.dict()
        progress_dict['completed'] = 0
        progress_dict['start_time'] = time.time()
        progress_dict['should_log'] = False
        progress_dict['last_log'] = None
        progress_dict['lock'] = manager.Lock()
        
        # Prepare arguments for each coordinate
        process_args = [
            (
                (y, x),
                weather_data,
                system,
                model_chain_kwargs,
                ptc,
                n_mods,
                progress_dict,
                idx,
                total_coords,
                compact_output,
            )
            for idx, (y, x) in enumerate(unique_coords, 1)
        ]
        
        coord_start_time = time.time()
        coord_subsets = []
        
        # Process coordinates in parallel
        if n_jobs == 1:
            # Sequential processing (useful for debugging or when only 1 coordinate)
            logger.debug("Using sequential processing (n_jobs=1)")
            for args in process_args:
                (y, x), subset = _process_single_coordinate(args)
                coord_subsets.append(subset)
                
                # Log progress
                idx = args[7]  # coord_index
                if idx % max(1, min(100, total_coords // 10)) == 0 or idx == total_coords:
                    elapsed_total = time.time() - coord_start_time
                    avg_time_per_coord = elapsed_total / idx
                    remaining_coords = total_coords - idx
                    eta = avg_time_per_coord * remaining_coords
                    logger.debug(
                        f"Processed coordinate {idx}/{total_coords} ({y:.2f}, {x:.2f}): "
                        f"Avg: {avg_time_per_coord:.2f}s/coord | "
                        f"ETA: {eta:.1f}s"
                    )
        else:
            # Parallel processing with progress tracking
            logger.debug(f"Using parallel processing with {n_jobs} workers")
            
            # Start a thread to monitor progress
            import threading
            stop_progress_thread = threading.Event()
            
            def progress_monitor():
                """Monitor progress and log updates"""
                last_logged = 0
                while not stop_progress_thread.is_set():
                    time.sleep(0.5)  # Check every 0.5 seconds
                    if progress_dict.get('should_log', False):
                        with progress_dict['lock']:
                            if progress_dict.get('should_log', False):
                                log_info = progress_dict.get('last_log')
                                if log_info and log_info['completed'] > last_logged:
                                    logger.debug(
                                        f"Processed coordinate {log_info['completed']}/{log_info['total']} "
                                        f"({log_info['coord'][0]:.2f}, {log_info['coord'][1]:.2f}): "
                                        f"Avg: {log_info['avg_time']:.2f}s/coord | "
                                        f"ETA: {log_info['eta']:.1f}s"
                                    )
                                    last_logged = log_info['completed']
                                progress_dict['should_log'] = False
            
            progress_thread = threading.Thread(target=progress_monitor, daemon=True)
            progress_thread.start()
            
            try:
                with Pool(processes=n_jobs) as pool:
                    results = pool.map(_process_single_coordinate, process_args)
                
                # Extract subsets from results
                coord_subsets = [subset for (y, x), subset in results]
                
            finally:
                stop_progress_thread.set()
                progress_thread.join(timeout=1.0)
        
        elapsed_total = time.time() - coord_start_time
        logger.info(
            f"Completed processing {total_coords} coordinate(s) in {elapsed_total:.2f}s "
            f"({elapsed_total/total_coords:.2f}s per coordinate on average)"
        )
        
        weather_data_final = pd.concat(coord_subsets).sort_index()

        out = xr.Dataset.from_dataframe(weather_data_final)
        desired_order = ("time", "x", "y")
        if all(d in out.dims for d in desired_order):
            out = out.transpose(*desired_order)
        return out
    
    def _prepare_dataset(self, source: xr.Dataset) -> xr.Dataset:
        """This will never be called, but must be implemented (abstract method)."""
        raise NotImplementedError("This model does not use _prepare_dataset")