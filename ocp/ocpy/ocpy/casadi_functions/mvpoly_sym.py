"""
mvpoly_sym.py

Symbolic multivariate polynomial basis and its Jacobian, built with
CasADi SX symbolics.

Equivalent to MATLAB's mvpoly_sym.m in the original codebase.
"""

from __future__ import annotations

import numpy as np
import casadi as ca
from itertools import combinations_with_replacement


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def mvpoly_sym(qs: ca.SX, order: int):
    """
    Build the monomial basis matrix and its partial-derivative matrix for
    a multivariate polynomial of given order.

    Parameters
    ----------
    qs : ca.SX, shape (1, n_dof)
        Symbolic row vector of joint coordinates.
    order : int
        Maximum total polynomial degree (inclusive).

    Returns
    -------
    mat : ca.SX, shape (n_points, n_coeff)
        Each row evaluates all monomials for one input point.
        Here n_points = qs.shape[0].
    diff_mat_q : ca.SX, shape (n_coeff, n_dof)
        diff_mat_q[:, j] holds the partial derivative of each monomial
        w.r.t. qs[j].
    """
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
    """
    Evaluate all monomials  x[0]^p0 * x[1]^p1 * ...  for each row of powers.

    x      : ca.SX, shape (1, n_dof)
    powers : np.ndarray, shape (n_coeff, n_dof)

    Returns ca.SX, shape (n_coeff, 1).
    """
    n_coeff, n_dof = powers.shape
    y = ca.SX.ones(n_coeff, 1)
    for col in range(n_dof):
        for row in range(n_coeff):
            p = int(powers[row, col])
            y[row] = y[row] * (x[0, col] ** p)
    return y


def _eval_der_monomial(x: ca.SX, powers: np.ndarray, idx: int) -> ca.SX:
    """
    Partial derivative of all monomials w.r.t. x[idx].

    d/dx_idx ( prod_j x_j^p_j ) = p_idx * x_idx^(p_idx-1) * prod_{j!=idx} x_j^p_j
    """
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
    """
    Return all monomial exponent tuples of total degree 0 … order.
    Shape: (n_coeff, n_dof).
    """
    rows = []
    for d in range(order + 1):
        rows.extend(_weak_compositions(d, n_dof))
    return np.array(rows, dtype=int)


def _weak_compositions(d: int, k: int) -> list[list[int]]:
    """
    All k-tuples of non-negative integers summing to d
    (equivalent to MATLAB's weak_compositions).
    """
    if k == 1:
        return [[d]]
    result = []
    for first in range(d + 1):
        for rest in _weak_compositions(d - first, k - 1):
            result.append([first] + rest)
    return result
