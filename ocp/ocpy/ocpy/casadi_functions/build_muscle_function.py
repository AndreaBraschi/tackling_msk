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
import sys, os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from casadi_functions.mvpoly_sym import mvpoly_sym


def build_muscle_function(path_muscle_poly: str | Path) -> ca.Function:
    """
    Build a CasADi Function for muscle geometry.

    Parameters
    ----------
    path_muscle_poly : str or Path
        Directory containing:
          - muscle_spanning_joints_info.mat
          - muscle_info_opt.mat

    Returns
    -------
    ca.Function
        'f_muscle' : (qin, qdotin) → (lMT, vMT, dM)
          qin     : (1, num_q)
          qdotin  : (1, num_q)
          lMT     : (num_muscles, 1)
          vMT     : (num_muscles, 1)
          dM      : (num_muscles, num_q)
    """
    path_muscle_poly = Path(path_muscle_poly)

    # ------------------------------------------------------------------ #
    # Load spanning-joint info and polynomial coefficients
    # ------------------------------------------------------------------ #
    spanning_mat  = scipy.io.loadmat(
        path_muscle_poly / "muscle_spanning_joints_info.mat"
    )
    muscle_info_mat = scipy.io.loadmat(
        path_muscle_poly / "muscle_info_opt.mat"
    )

    # The first (non-meta) field holds the spanning-joint matrix
    field = next(k for k in spanning_mat if not k.startswith("_"))
    spanning = spanning_mat[field]                  # (num_muscles, num_q)

    muscle_info = muscle_info_mat["MuscleInfo"]     # structured array

    num_q, num_muscles = spanning.shape[1], spanning.shape[0]

    # ------------------------------------------------------------------ #
    # Symbolic design variables
    # ------------------------------------------------------------------ #
    qin    = ca.SX.sym("qin",    1, num_q)
    qdotin = ca.SX.sym("qdotin", 1, num_q)

    lMT = ca.SX(num_muscles, 1)
    vMT = ca.SX(num_muscles, 1)
    dM  = ca.SX(num_muscles, num_q)

    for i in range(num_muscles):
        # DoF indices this muscle spans (0-based)
        dof_indices = np.where(spanning[i, :] == 1)[0]

        # Polynomial order for this muscle
        order = int(muscle_info["muscle"][0, 0]["order"][0, i])

        # Coefficients stored in a cell array; pick the one for 'order'
        # In scipy, MATLAB cell arrays come out as object arrays
        coeff_cell = muscle_info["muscle"][0, 0]["coeff"][0, i]
        # coeff_cell is shape (1, order); pick last (index = order-1)
        coeff = np.squeeze(coeff_cell[0, order - 1])   # (n_coeff,)

        # Symbolic polynomial basis
        q_subset = qin[0, dof_indices]                 # (1, n_spanning)
        mat, diff_mat_q = mvpoly_sym(q_subset, order)  # (1,n_coeff), (n_coeff,n_sp)

        # Non-zero coefficient indices
        nz = np.where(coeff != 0)[0]
        coeff_nz = coeff[nz]

        # Muscle-tendon length
        lMT[i] = mat[0, nz] @ coeff_nz

        vMT[i] = 0
        dM[i, :] = 0

        for dof_nr, glob_idx in enumerate(dof_indices):
            # Moment arm: –d(lMT)/d(q)
            ma = -(diff_mat_q[nz, dof_nr].T @ coeff_nz)
            dM[i, glob_idx] = ma

            # Rate of change of lMT
            vMT[i] = vMT[i] + (-dM[i, glob_idx] * qdotin[0, glob_idx])

    return ca.Function("f_muscle", [qin, qdotin], [lMT, vMT, dM])
