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


import json

import numpy as np
import pandas as pd
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
        maybe_progressbar = lambda x: x  # noqa: E731

    return maybe_progressbar


def dummy_njit(f=None, *args, **kwargs):
    """Dummy decorator for numba.njit. Handles the case when numba is not installed.

    Args:
        f (function): Function to be decorated. If None, returns identity.
    """

    def decorator(func):
        return func

    if callable(f):
        return f

    return decorator


def get_daterange(years: slice, months: slice):
    """Get the date range covering the entire years and months range.

    Args:
        years (slice): The years range.
        months (slice): The months range.

    Returns:
        pd.
    """

    assert years.start <= years.stop, "Start year must be less than stop year."
    assert months.start <= months.stop, "Start month must be less than stop month."

    return pd.date_range(
        start=pd.Timestamp(f"{years.start}-{months.start}-1"),
        end=pd.Timestamp(
            f"{years.stop}-{months.stop}-{pd.Timestamp(f'{years.stop}-{months.stop}-1').days_in_month}"
        ),
        freq="d",
    )


class NpEncoder(json.JSONEncoder):
    """Custom JSON encoder for NumPy data types, as the default JSON encoder does not handle them."""

    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.floating):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return super(NpEncoder, self).default(obj)
