"""
build_muscle_function.py

Builds a CasADi Function that computes, symbolically, the muscle-tendon
lengths (lMT), velocities (vMT), and moment arms (dM) from joint
coordinates and velocities using pre-fitted polynomial coefficients.
"""

import scipy.io
from pathlib import Path
import casadi as ca
import numpy as np


from .mvpoly_sym import mvpoly_sym


def build_muscle_function(muscle_pkl) -> ca.Function:

    # which joints are spanned by which muscle
    muscle_spanning = muscle_pkl['spanning']
    # where polynomial infos are stored
    muscle_info = muscle_pkl['muscle_info']

    num_q, num_muscles = muscle_spanning.shape[1], muscle_spanning.shape[0]

    # Symbolic variables for subset of q's and q dot needed to compute muscle lengths and moment arms
    qin    = ca.SX.sym("qin", 1, num_q)
    qdotin = ca.SX.sym("qdotin", 1, num_q)

    lMT = ca.SX(num_muscles, 1)             # muscle-tendon length
    vMT = ca.SX(num_muscles, 1)             # rate of change of muscle-tendon length
    dM  = ca.SX(num_muscles, num_q)         # muscle-tendon moment arm

    for i, muscle in enumerate(muscle_info):

        # find which DoFs the current muscle spans
        dof_indices = np.where(muscle_spanning[i, :] == 1)[0].tolist()

        # Polynomial order for this muscle
        order = muscle_info[muscle]['order']

        # poly coefficients
        coeff = muscle_info[muscle]['coeff_selected']

        # retrieve only the polynomial that were selected after the reduction (those will be flagged by a 0 in the array
        # where the coefficients are stored).
        indices_nonzero = np.where(coeff != 0)[0]
        coeff_nonzero = coeff[indices_nonzero]

        # compute muscle length and moment arm, symbolically, given the set of
        # q's that the current muscle spans.
        mat, diff_mat_q = mvpoly_sym(qin[:, dof_indices], order)  # (1,n_coeff), (n_coeff,n_sp)


        # compute muscle-tendon length only using the relevant coefficients
        lMT[i] = mat[0, indices_nonzero] @ coeff_nonzero

        vMT[i] = 0
        dM[i, :] = 0

        for qi_local, qi_global in enumerate(dof_indices):
            # Moment arm: –d(lMT)/d(q)
            ma = -(diff_mat_q[indices_nonzero, qi_local].T @ coeff_nonzero)
            dM[i, qi_global] = ma


            vMT[i] = vMT[i] + (-dM[i, qi_global] * qdotin[0, qi_global])

    return ca.Function("f_muscle", [qin, qdotin], [lMT, vMT, dM])
