from __future__ import annotations

import numpy as np
import casadi as ca


def mvpoly_sym(qs: ca.SX, order: int):

    n_dof    = qs.shape[1]
    n_points = qs.shape[0]

    exponents = _all_monomials_upto(order, n_dof)   # (n_coeff, n_dof)
    n_coeff   = exponents.shape[0]

    mat = ca.SX(n_points, n_coeff)
    for i in range(n_points):
        mat[i, :] = _eval_monomial(qs[i, :], exponents).T

    diff_mat_q = ca.SX(n_coeff, n_dof)
    for j in range(n_dof):
        for i in range(n_points):
            diff_mat_q[:, j] = _eval_der_monomial(qs[i, :], exponents, j)

    return mat, diff_mat_q


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _eval_monomial(x: ca.SX, powers: np.ndarray) -> ca.SX:

    n_coeff, n_dof = powers.shape
    y = ca.SX.ones(n_coeff, 1)
    for col in range(n_dof):
        for row in range(n_coeff):
            p = int(powers[row, col])
            y[row] = y[row] * (x[0, col] ** p)
    return y


def _eval_der_monomial(x: ca.SX, powers: np.ndarray, idx: int) -> ca.SX:

    n_coeff, n_dof = powers.shape
    coeff   = powers[:, idx].copy().astype(float)   # p_idx values
    new_pow = powers.copy().astype(float)
    new_pow[:, idx] -= 1
    new_pow[new_pow < 0] = 0                         # clamp (avoids 0^-1)

    y = ca.SX(n_coeff, 1)
    for row in range(n_coeff):
        val = coeff[row]
        for col in range(n_dof):
            val = val * (x[0, col] ** int(new_pow[row, col]))
        y[row] = val
    return y


def _all_monomials_upto(order: int, n_dof: int) -> np.ndarray:

    rows = []
    for d in range(order + 1):
        rows.extend(_weak_compositions(d, n_dof))
    return np.array(rows, dtype=int)


def _weak_compositions(d: int, k: int) -> list[list[int]]:

    if k == 1:
        return [[d]]
    result = []
    for first in range(d + 1):
        for rest in _weak_compositions(d - first, k - 1):
            result.append([first] + rest)
    return result
