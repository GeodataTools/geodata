# Copyright 2023 Michael Davidson (UCSD), Xiqiang Liu (UCSD)

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


# pylint: disable=unused-argument
def dummy_njit(f=None, *args, **kwargs):  #  pylint: disable=keyword-arg-before-vararg
    """Dummy decorator for numba.njit. Handles the case when numba is not installed.

    Args:
        f (function): Function to be decorated. If None, returns identity.
    """

    def decorator(func):
        return func

    if callable(f):
        return f

    return decorator


# pylint: enable=unused-argument
