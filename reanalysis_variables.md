Variables
-


## Variables used in calculations

**Variables used for wind calculation (MERRA/MERRA2)**
ustar
z0m
disph
rhoa
ulml
vlml
tstar
hlml
tlml
pblh
hflux


**Variables used for solar calculation (ERA5)**
(shortName) | (name)                        | (paramId)
z           | Geopotential (CDS: Orography) | 129
tisr        | TOA incident solar radiation                | 212
ssrd        | Surface Solar Rad Downwards                 | 169
ssr         | Surface net Solar Radiation                 | 176
fdir        | Total sky direct solar radiation at surface | 228021
ro          | Runoff                                      | 205
2t          | 2 metre temperature                         | 167
sp          | Surface pressure                            | 134
stl4        | Soil temperature level 4                    | 236
fsr         | Forecast surface roughness                   | 244


## Files

**MERRA-2**

tavg1_2d_flx_Nx (M2T1NXFLX): Surface Flux Diagnostics
- all the usual
- Granule Size: ~379 MB


tavg1_2d_rad_Nx (M2T1NXRAD): Radiation Diagnostics
- albedo
- toa incoming shortwave flux
- ...
- Granule Size: ~209 MB


**MERRA**

tavg1_2d_flx_Nx
- Granule size: ~ 277 MB



**Misc notes on computation**

Deleting variables: http://www.cgd.ucar.edu/cms/nco
https://yidongwonyi.wordpress.com/linux-data-handling-netcdf-nc/nco-extract-variable-delete-variabledimension/
https://stackoverflow.com/questions/20215529/delete-a-dimension-in-a-netcdf-file

Chunking and deflating helps with I/O: http://wiki.seas.harvard.edu/geos-chem/index.php/Working_with_netCDF_data_files#Chunking_and_deflating_a_netCDF_file_to_improve_I.2FO
