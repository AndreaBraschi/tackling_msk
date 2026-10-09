"""
eval_spline.py

Evaluates a cubic spline (and optionally its 1st and 2nd derivatives)
at a set of points x, given a PPoly spline object (scipy equivalent of
MATLAB's ppval/spline).

NOTE: This function assumes a cubic spline (degree 3). Derivatives are
computed analytically from the spline coefficients.

Analytical expressions:
  1st derivative: y_dot  =  3*a1*dx^2 + 2*a2*dx + a3
  2nd derivative: y_ddot =  6*a1*dx   + 2*a2
where dx is the offset within each spline segment.
"""

import numpy as np
from scipy.interpolate import CubicSpline
from typing import Tuple

def eval_spline(cs: CubicSpline, x: np.ndarray, flag: int) -> Tuple[np.ndarray, ...]:
    """
    Evaluate a cubic spline and, optionally, its 1st and 2nd derivatives.

    Parameters
    ----------
    cs : scipy.interpolate.PPoly
        Piecewise-polynomial cubic spline (e.g. from scipy.interpolate.CubicSpline).
    x : np.ndarray
        1-D array of evaluation points.
    flag : int
        0 → position only
        1 → position + velocity
        2 → position + velocity + acceleration

    Returns
    -------
    dict with keys 'pos', and optionally 'vel', 'acc'.
    """

    pos = cs(x)

    if flag >= 1:
        vel = _cubic_spline_vel(cs, x)
    if flag >= 2:
        acc = _cubic_spline_acc(cs, x)

    return pos, vel, acc


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _cubic_spline_vel(cs: CubicSpline, x: np.ndarray) -> np.ndarray:
    """
    Velocity from a cubic spline using analytical differentiation.

    MATLAB stores coefficients per breakpoint interval as:
        p(t) = a1*(t-t_k)^3 + a2*(t-t_k)^2 + a3*(t-t_k) + a4
    so the derivative is:
        p'(t) = 3*a1*(t-t_k)^2 + 2*a2*(t-t_k) + a3

    At the *left* end of each interval dt=0, so the velocity = a3 (index 2).
    At the very last point we compute it explicitly for dx = x[-1] - x[-2].
    """
    # cs.c has shape [order, n_intervals, num_dimensions]; cubic → order = 4
    # coefficients: cs.c[0]=a1, cs.c[1]=a2, cs.c[2]=a3, cs.c[3]=a4

    T = cs.c.shape[-2] + 1    # we add 1 as from cs.c we know the number of intervals, but we need the number of knots
    dims = cs.c.shape[-1]

    a1 = cs.c[0]   # cubic coefficient
    a2 = cs.c[1]   # quadratic coefficient
    a3 = cs.c[2]   # linear coefficient (= derivative at left knot)

    vel = np.zeros((T, dims))
    # Interior / left-end points: at dt=0 the derivative is simply a3
    vel[:-1] = a3
    # Last point: evaluate analytically at dt = x[-1] - x[-2]
    h = x[-1] - x[-2]
    vel[-1] = 3 * a1[-1] * h**2 + 2 * a2[-1] * h + a3[-1]

    return vel


def _cubic_spline_acc(cs: CubicSpline, x: np.ndarray) -> np.ndarray:
    """
    Acceleration from a cubic spline using analytical differentiation.

    p''(t) = 6*a1*(t-t_k) + 2*a2

    At the left end of each interval dt=0, so acceleration = 2*a2.
    At the very last point we compute it explicitly.
    """
    T = cs.c.shape[-2] + 1
    dims = cs.c.shape[-1]

    a1 = cs.c[0]
    a2 = cs.c[1]

    acc = np.zeros((T, dims))
    acc[:-1] = 2 * a2
    h = x[-1] - x[-2]
    acc[-1] = 6 * a1[-1] * h + 2 * a2[-1]

    return acc
