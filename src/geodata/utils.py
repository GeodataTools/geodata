# Copyright 2023, 2025 Michael Davidson (UCSD), Xiqiang Liu (UCSD)

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


import hashlib
import json
from pathlib import Path

import dask.array as da
import numpy as np
import pandas as pd
import xarray as xr


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
        pd.DatetimeIndex: The date range.
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


def ensure_slice(obj: slice | list):
    """Ensure that the input is a slice object. If the input is a list, convert it to a slice object.

    Args:
        obj (slice | list): The input object.

    Returns:
        slice: The converted slice object.
    """
    if isinstance(obj, list) and (len(obj) == 2 or len(obj) == 3):
        return slice(*obj)
    elif isinstance(obj, slice):
        return obj
    else:
        raise TypeError("Input must be a slice or a list.")


def check_hash(file: Path, saved_hash: str | None = None) -> tuple[bool, str]:
    """Check if the hash of a file matches the given hash.

    Args:
        file (Path): The path to the file.
        saved_hash (str | None): The hash to compare against. If None, the function will compute the hash of the file.
            If provided, the function will compare the computed hash with this value.
            If the hashes match, the function will return True and the computed hash.
            If the hashes do not match, the function will return False and the computed hash.

    Returns:
        tuple[bool, str]: A tuple containing a boolean indicating if the hash matches and the computed hash.
    """

    if not file.exists():
        return False, ""

    with open(file, "rb") as f:
        computed_hash = hashlib.sha256(f.read()).hexdigest()

    if saved_hash is None:
        return True, computed_hash
    else:
        return computed_hash == saved_hash, computed_hash


def rechunk_dataset(
    data: xr.DataArray | xr.Dataset,
    target_chunk_bytes: int = 20 * 1024**2,
    force_full_chunk_dims: list[str] | None = None,
):
    """
    Rechunk xarray DataArray or Dataset to maximize chunk size
    under a memory limit, while forcing certain dimensions to be unchunked
    (i.e., use only one chunk across that dimension).

    Parameters:
        data: xr.DataArray or xr.Dataset
        target_chunk_bytes: maximum memory per chunk (in bytes)
        force_full_chunk_dims: list of dimension names to not chunk (single chunk along that dim)
    """
    if force_full_chunk_dims is None:
        force_full_chunk_dims = []

    if isinstance(data, xr.Dataset):
        vars_to_chunk = {
            name: rechunk_dataset(var, target_chunk_bytes, force_full_chunk_dims)
            for name, var in data.data_vars.items()
        }
        return data.assign(vars_to_chunk)

    if not isinstance(data.data, da.Array):
        raise ValueError("Data must be a Dask-backed xarray object")

    shape = data.shape
    dims = data.dims
    itemsize = data.dtype.itemsize

    # Start with full dims
    chunk_shape = list(shape)
    dim_to_index = {dim: i for i, dim in enumerate(dims)}

    # Force full chunks on specified dims
    for dim in force_full_chunk_dims:
        if dim in dim_to_index:
            chunk_shape[dim_to_index[dim]] = shape[dim_to_index[dim]]

    # Reduce non-fixed dims to fit memory budget
    while True:
        est_bytes = np.prod(chunk_shape) * itemsize
        if est_bytes <= target_chunk_bytes:
            break

        # Pick largest non-fixed dimension to halve
        candidates = [
            (i, size)
            for i, size in enumerate(chunk_shape)
            if dims[i] not in force_full_chunk_dims and size > 1
        ]
        if not candidates:
            break  # Can't reduce further

        i, _ = max(candidates, key=lambda x: x[1])
        chunk_shape[i] = max(1, chunk_shape[i] // 2)

    chunk_dict = dict(zip(dims, chunk_shape))
    return data.chunk(chunk_dict)
