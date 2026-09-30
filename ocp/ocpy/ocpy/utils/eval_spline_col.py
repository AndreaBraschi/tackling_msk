"""
eval_spline_col.py

Like eval_spline.py but designed for collocation points where we already
know which spline segment (index) each point belongs to and the local
offset (dt) within that segment.

This avoids a redundant binary search inside ppval / PPoly.__call__.

Analytical expressions (cubic spline):
  1st derivative: y_dot  = 3*a1*dt^2 + 2*a2*dt + a3
  2nd derivative: y_ddot = 6*a1*dt   + 2*a2
"""

import numpy as np
from scipy.interpolate import CubicSpline
from typing import Tuple

def eval_spline_col(cs: CubicSpline, x: np.ndarray, indices: np.ndarray, dt: np.ndarray) -> Tuple[np.ndarray, ...]:
    """
    Evaluate a cubic spline at collocation points.

    Parameters
    ----------
    cs : scipy.interpolate.PPoly
        Piecewise-polynomial cubic spline.
    x : np.ndarray
        Evaluation points (used for the position via cs(x)).
    indices : np.ndarray
        Integer array – for each point in x, the index of the spline
        interval (0-based) it belongs to.
    dt : np.ndarray
        Local time offset within the interval: dt = x - breakpoints[indices].
    flag : int
        0 → position only; 1 → + velocity; 2 → + velocity + acceleration.

    Returns
    -------
    dict with keys 'pos', and optionally 'vel', 'acc'.
    """
    pos = cs(x)
    vel = _cubic_spline_vel_col(cs, indices, dt)
    acc = _cubic_spline_acc_col(cs, indices, dt)

    return pos, vel, acc


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _cubic_spline_vel_col(cs: CubicSpline, indices: np.ndarray, dt: np.ndarray) -> np.ndarray:
    """
    Velocity using pre-computed segment indices and local offsets.

    p'(t) = 3*a1*dt^2 + 2*a2*dt + a3
    """
    a1 = cs.c[0]   # cubic coefficient
    a2 = cs.c[1]   # quadratic coefficient
    a3 = cs.c[2]   # linear coefficient (= derivative at left knot)

    vel = 3 * a1[indices] * dt**2 + 2 * a2[indices] * dt + a3[indices]
    return vel


def _cubic_spline_acc_col(cs: CubicSpline, indices: np.ndarray, dt: np.ndarray) -> np.ndarray:
    """
    Acceleration using pre-computed segment indices and local offsets.

    p''(t) = 6*a1*dt + 2*a2
    """
    a1 = cs.c[0]
    a2 = cs.c[1]

    acc = 6 * a1[indices] * dt + 2 * a2[indices]
    return acc
