"""
hill_equilibrium.py

Builds a CasADi Function that encodes the Hill-type muscle force
equilibrium (De Groote et al. 2016).
"""

import casadi as ca
import numpy as np
from ..dynamics_constraints.muscle_force_equilibrium import force_equilibrium


def generate_hill_equilibrium_func(num_muscles: int, MT_parameters: np.ndarray) -> ca.Function:


    # Symbolic inputs
    fT  = ca.SX.sym("fT",  num_muscles)        # normalised tendon forces
    a_m  = ca.SX.sym("am", num_muscles)        # muscle activations
    dfT = ca.SX.sym("dfT", num_muscles)        # normalised tendo force time derivatives
    lMT = ca.SX.sym("lMT", num_muscles)        # MTU lengths
    vMT = ca.SX.sym("vMT", num_muscles)        # MTU length rate of change

    # Symbolic outputs  (pre-allocate)
    hill_diff = ca.SX(num_muscles, 1)

    for m in range(num_muscles):
        # Slice single-muscle parameters: shape (5, 1)
        params_m = MT_parameters[:, m : m + 1]
        hill_diff[m] = force_equilibrium(
            a_m[m], fT[m], dfT[m], lMT[m], vMT[m], params_m)

    return ca.Function(
        "hill_equilibrium",
        [a_m, fT, dfT, lMT, vMT],
        [hill_diff]
    )
