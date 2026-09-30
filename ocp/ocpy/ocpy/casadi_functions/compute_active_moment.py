"""
compute_active_moment.py

Builds a CasADi Function that computes the active joint moment as the
dot product of moment arms (r) and muscle forces (force):
    M = sum(r * force)
"""

import casadi as ca


def compute_active_moment(D: int) -> ca.Function:
    """
    Create a CasADi Function for the active joint moment.

    Parameters
    ----------
    D : int
        Number of muscles contributing to the moment.

    Returns
    -------
    ca.Function
        CasADi Function 'M' that maps (r, force) → scalar M.
    """
    r     = ca.SX.sym("r",     D)
    force = ca.SX.sym("force", D)

    M = ca.sum1(r * force)

    return ca.Function("M", [r, force], [M])
