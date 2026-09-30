"""
read_mot_grf.py

Reads a .mot file containing Ground Reaction Forces (GRFs) and returns
a dictionary equivalent to MATLAB's struct output.

A Butterworth low-pass filter is optionally applied when a cutoff
frequency (Hz) is provided.
"""

import numpy as np
from scipy.signal import butter, filtfilt
from pathlib import Path


def read_mot_grf(filepath: str | Path, cutoff: float | None = None) -> dict:
    """
    Read a .mot GRF file.

    Parameters
    ----------
    filepath : str or Path
        Full path to the .mot GRF file.
    cutoff : float or None
        Low-pass filter cutoff frequency in Hz. If None (default), no
        filtering is applied.

    Returns
    -------
    dict with keys:
        'data'      – (N, M) float array; column 0 is time.
        'colheaders'– list of column-header strings.
        'time'      – (N,) float array (= data[:, 0]).
    """
    filepath = Path(filepath)

    col_headers, data = _parse_mot_file(filepath)

    grf: dict = {
        "colheaders": col_headers,
        "data": data,
        "time": data[:, 0],
    }

    if cutoff is not None:
        order = 2
        fs = 1.0 / float(np.mean(np.diff(data[:, 0])))
        b, a = butter(order // 2, cutoff / (0.5 * fs), btype="low")
        grf_filt = filtfilt(b, a, data[:, 1:], axis=0)
        grf["data"] = np.hstack([data[:, :1], grf_filt])

    return grf


# ---------------------------------------------------------------------------
# Internal helper
# ---------------------------------------------------------------------------

def _parse_mot_file(filepath: Path):
    """
    Parse an OpenSim .mot file, returning (col_headers, data_array).

    .mot files contain a text header terminated by 'endheader', followed
    by a tab-separated numerical block whose first row is the column names.
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
