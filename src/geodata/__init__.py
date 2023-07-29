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
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

from ._version import __version__
from .cutout import Cutout
from .dataset import Dataset
from .mask import Mask
from .plot import *  # noqa: F403
from .model import *  # noqa: F403

__author__ = "Michael Davidson (UCSD), William Honaker"
__copyright__ = "GNU GPL 3 license"


__all__ = ["Cutout", "Dataset", "Mask", "__version__"]
