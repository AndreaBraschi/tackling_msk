"""
read_mot_qs.py

Reads a .mot IK result file, converts rotational DoFs from degrees to
radians, and optionally applies a Butterworth low-pass filter.

Returns a dictionary equivalent to MATLAB's Qs struct.
"""

import numpy as np
from scipy.signal import butter, filtfilt
from pathlib import Path


# Translational DoF suffixes that should NOT be converted to radians.
_TRANSLATIONAL_TAGS = ("tx", "ty", "tz")


def read_mot_qs(filepath: str | Path, cutoff: float | None = None) -> dict:
    """
    Read a .mot IK file.

    Parameters
    ----------
    filepath : str or Path
        Full path to the .mot file.
    cutoff : float or None
        Low-pass filter cutoff frequency in Hz. If None (default), no
        filtering is applied.

    Returns
    -------
    dict with keys:
        'data'         – (N, M) float array; column 0 is time (seconds),
                         rotational columns in radians.
        'colheaders'   – list of column-header strings.
        'time'         – (N,) float array (= data[:, 0]).
        'allfilt'      – filtered copy of 'data' (or unfiltered copy when
                         cutoff is None).
    """
    filepath = Path(filepath)

    col_headers, data = _parse_mot_file(filepath)

    # Build a boolean mask: True for columns that are translational (keep degrees)
    mask = np.array(
        [any(tag in h for tag in _TRANSLATIONAL_TAGS) for h in col_headers],
        dtype=bool,
    )

    # Convert rotational columns from degrees → radians (skip time col at index 0)
    data[:, ~mask] = data[:, ~mask] * (np.pi / 180.0)

    qs: dict = {
        "colheaders": col_headers,
        "data": data,
        "time": data[:, 0],
    }

    if cutoff is not None:
        order = 2
        fs = 1.0 / float(np.mean(np.diff(data[:, 0])))
        b, a = butter(order // 2, cutoff / (0.5 * fs), btype="low")
        allfilt = data.copy()
        allfilt[:, 1:] = filtfilt(b, a, data[:, 1:], axis=0)
        qs["allfilt"] = allfilt
    else:
        qs["allfilt"] = data.copy()

    return qs


# ---------------------------------------------------------------------------
# Internal helper
# ---------------------------------------------------------------------------

def _parse_mot_file(filepath: Path):
    """
    Parse an OpenSim .mot file → (col_headers, data_array).
    """
    col_headers = []
    data_rows = []
    in_data = False

    with open(filepath, "r") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not in_data:
                if line.strip().lower() == "endheader":
                    in_data = True
            else:
                if not col_headers:
                    col_headers = line.split()
                else:
                    tokens = line.split()
                    if tokens:
                        data_rows.append([float(t) for t in tokens])

    return col_headers, np.array(data_rows, dtype=float)
