# Copyright 2024 Michael Davidson (UCSD), Xiqiang Liu (UCSD)

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

from herbie import Herbie

from ._base import HRRRBaseDataset


class HRRRWindDataset(HRRRBaseDataset):
    """
    HRRRWindDataset is a class that encaps a dataset from the HRRR
    dataset. It provides a streamlined workflow for downloading, preprocessing,
    and storing of these datasets.

    The HRRR dataset is a high-resolution weather forecast model that provides
    hourly data for the United States. This dataset is useful for a variety of
    applications, including renewable energy forecasting, weather prediction,
    and climate research.

    This class provides a simple interface for downloading and processing the
    HRRR dataset. It allows users to specify the years and months of interest,
    as well as the variables they wish to download.
    """
