# Copyright 2025 Michael Davidson (UCSD), Xiqiang Liu (UCSD)

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


from .._base import MERRA2BaseDataset


class MERRA2SLVRadiationHourlyDataset(MERRA2BaseDataset):
    """MERRA2SLVRadiationHourlyDataset is a class that encaps a dataset from the MERRA2 reanalysis
    dataset. It provides a streamlined workflow for downloading, preprocessing,
    and storing of these datasets.
    """

    weather_config = "slv_radiation_hourly"

    variables = ["albedo", "swgdn", "swtdn", "t2m"]

    url_template = (
        "https://goldsmr4.gesdisc.eosdis.nasa.gov/data/MERRA2/M2T1NXSLV.5.12.4/{year}/{month:0>2}/MERRA2_{spinup}.tavg1_2d_slv_Nx.{year}{month:0>2}{day:0>2}.nc4",
        "https://goldsmr4.gesdisc.eosdis.nasa.gov/data/MERRA2/M2T1NXRAD.5.12.4/{year}/{month:0>2}/MERRA2_{spinup}.tavg1_2d_rad_Nx.{year}{month:0>2}{day:0>2}.nc4",
    )
