# MERRA2 Configs and Outputs

A list of currently possible **geodata** configs and outputs for MERRA2 data from NASA's GES DISC.

## Dataset Configurations

**geodata** is currently optimized to work with the following MERRA2 datasets:

* `surface_flux_hourly`: [MERRA2 hourly, single-level surface flux diagnostics](https://disc.gsfc.nasa.gov/datasets/M2T1NXFLX_5.12.4/summary)
* `surface_flux_monthly`: [MERRA2 monthly mean, single-level surface flux diagnostics](https://disc.gsfc.nasa.gov/datasets/M2TMNXFLX_5.12.4/summary)
* `surface_flux_dailymeans`: [MERRA2 daily mean, single-level diagnostics](https://disc.gsfc.nasa.gov/datasets/M2SDNXSLV_5.12.4/summary)
* `slv_radiation_hourly`: 
    - [MERRA2 hourly, single-level diagnostics](https://disc.gsfc.nasa.gov/datasets/M2SDNXSLV_5.12.4/summary)
    - [MERRA2 hourly, single-level radiation diagnostics](https://disc.gsfc.nasa.gov/datasets/M2T1NXRAD_5.12.4/summary)
* `slv_radiation_monthly`:
    - [MERRA2 monthly mean, single-level diagnostics](https://disc.gsfc.nasa.gov/datasets/M2TMNXSLV_5.12.4/summary)
    - [MERRA2 monthly mean, single-level radiation diagnostics](https://disc.gsfc.nasa.gov/datasets/M2TMNXRAD_5.12.4/summary)


## Variable Configurations

MERRA2 surface flux data currently supports the following variable configs:

**Wind** (`surface_flux_hourly`, `surface_flux_monthly`):
* ustar: surface_velocity_scale (m s-1)
* z0m: surface_roughness (m)
* disph: zero_plane_displacement_height (m)
* rhoa: air_density_at_surface (kg m-3)
* ulml: surface_eastward_wind (m s-1)
* vlml: surface_northward_wind (m s-1)
* tstar: surface_temperature_scale (K)
* hlml: surface_layer_height (m)
* tlml: surface_air_temperature (K)
* pblh: planetary_boundary_layer_height (m)
* hflux: sensible_heat_flux_from_turbulence (W m-2)
* eflux: total_latent_energy_flux (W m-2)

**Solar** (`slv_radiation_hourly`, `slv_radiation_monthly`)
* albedo: surface albedo 
* swgdn: surface incoming shortwave flux  (W m-2)
* swtdn: toa incoming shortwave flux  (W m-2)
* t2m: 2-meter air temperature (K)

**Temperature** (`surface_flux_dailymeans`)
* hournorain: time-during an hour with no precipitation (s)
* tprecmax: Maximum precipitation rate during the period (kg m-2 s-1)
* t2max: max 2-meter air temperature (K)
* t2mmean: mean 2-meter air temperature (K)
* t2mmin: min 2-meter air temperature (K)

## Outputs

MERRA2 surface flux data currently supports the following outputs:

**Wind** (`surface_flux_hourly`, `surface_flux_monthly`):
* [Wind generation time-series]
* [Wind speed time-series]
* [Wind power density time-series]

**Solar** (`slv_radiation_hourly`, `slv_radiation_monthly`)
* [pv generation time-series]

**Temperature** (`surface_flux_dailymeans`)
* [Celsius Temperature]
