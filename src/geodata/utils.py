## Copyright 2016-2017 Gorm Andresen (Aarhus University), Jonas Hoersch (FIAS), Tom Brown (FIAS)

## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License as
## published by the Free Software Foundation; either version 3 of the
## License, or (at your option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.

## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.


"""
GEODATA

Geospatial Data Collection and "Pre-Analysis" Tools
"""

import progressbar as pgb


def make_optional_progressbar(show, prefix, max_value):
    if show:
        widgets = [
            pgb.widgets.Percentage(),
            " ",
            pgb.widgets.SimpleProgress(),
            " ",
            pgb.widgets.Bar(),
            " ",
            pgb.widgets.Timer(),
            " ",
            pgb.widgets.ETA(),
        ]
        if not prefix.endswith(": "):
            prefix = prefix.strip() + ": "
        maybe_progressbar = pgb.ProgressBar(
            prefix=prefix, widgets=widgets, max_value=max_value
        )
    else:
        maybe_progressbar = lambda x: x  # pylint: disable=unnecessary-lambda-assignment

    return maybe_progressbar


def dummy_njit(f=None, *args, **kwargs):  #  pylint: disable=keyword-arg-before-vararg
    """Dummy decorator for numba.njit. Handles the case when numba is not installed.

    Args:
        f (function): Function to be decorated. If None, returns a an identity decorator.
    """

    def decorator(func):
        return func

    if callable(f):
        return f
    else:
        return decorator
