# MERRA2 Configs and Outputs

## Variable Configurations

MERRA2 surface flux data currently supports the following variable configs:

**Wind** (`surface_flux`):
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

**Solar** (`TBD`)

## Outputs

MERRA2 surface flux data currently supports the following outputs:

**Wind**
* [Wind speed for turbine to wind energy generation]
* Wind speeds at given height
* Wind power density at given height (WPD = 0.5 * Density * Windspd^3)
* Wind generation time-series
* Wind speed time-series
* Wind power density time-series