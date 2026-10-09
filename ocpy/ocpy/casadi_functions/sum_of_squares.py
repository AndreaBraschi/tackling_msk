"""
sum_of_squares.py

Builds a CasADi Function that computes the sum of squares of a symbolic
vector of dimension D.
"""

import casadi as ca


def sum_of_squares(var_name: str, D: int) -> ca.Function:
    """
    Create a CasADi Function for  J = sum(var^2).

    Parameters
    ----------
    var_name : str
        Name given to the symbolic variable (used in the expression graph).
    D : int
        Dimensionality (length) of the vector.

    Returns
    -------
    ca.Function
        CasADi Function named 'J_<var_name>' that maps var → scalar J.
    """
    var = ca.SX.sym(var_name, D)
    J = ca.sum1(var**2)
    return ca.Function(f"J_{var_name}", [var], [J])
