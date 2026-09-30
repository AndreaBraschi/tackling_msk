"""
muscle_force_equilibrium.py

Hill-type muscle model: De Groote et al. (2016) formulation.

Original author: Antoine Falisse (12/19/2018)
Reference: De Groote et al., Ann Biomed Eng (2016)
           DOI: 10.1007/s10439-016-1591-9

Works with NumPy arrays and CasADi symbolics alike (all operations are
element-wise and use only arithmetic/exp/log/sqrt, which CasADi overloads).
"""

import numpy as np


def force_equilibrium(a, fse, dfse, lMT, vMT, params, Fvparam, Fpparam, Faparam):
    """
    Compute the Hill-equilibrium and related muscle quantities.

    Parameters
    ----------
    a      : activation  (num_muscles,)
    fse    : normalised tendon force  FTtilde  (num_muscles,)
    dfse   : time derivative of normalised tendon force  (num_muscles,)
    lMT    : muscle-tendon length  (num_muscles,)
    vMT    : muscle-tendon velocity  (num_muscles,)
    params : (5, num_muscles) float array
                row 0: FMo   – max isometric force
                row 1: lMo   – optimal fibre length
                row 2: lTs   – tendon slack length
                row 3: alphao– pennation angle at optimal fibre length
                row 4: vMmax – max contraction velocity (= vMaxrel * lMo)
    Fvparam : (4,) force-velocity curve parameters
    Fpparam : (2,) passive force-length curve parameters
    Faparam : (8,) active force-length curve parameters

    Returns
    -------
    err      : Hill equilibrium error  (Fce + Fpe)*cos_alpha - fse
    FT       : tendon force            = fse * FMo
    Fce      : contractile element force
    Fiso     : normalised active force-length value  FMltilde
    vMmax    : max contraction velocity
    Fpetilde : normalised passive force
    lMtilde  : normalised fibre length
    """
    # Unpack muscle-tendon parameters (broadcast over muscles)
    FMo   = params[0, :]
    lMo   = params[1, :]
    lTs   = params[2, :]
    alphao = params[3, :]
    vMmax  = params[4, :]

    Atendonsc = 35.0
    Atendon = Atendonsc  # scalar; broadcast will handle arrays

    # ------------------------------------------------------------------ #
    # Inverse tendon force-length: recover normalised tendon length lTtilde
    # ------------------------------------------------------------------ #
    lTtilde = _log(5.0 * (fse + 0.25)) / Atendon + 0.995

    # ------------------------------------------------------------------ #
    # Geometric relationship: fibre length
    # ------------------------------------------------------------------ #
    lM      = _sqrt((lMo * _sin(alphao))**2 + (lMT - lTs * lTtilde)**2)
    lMtilde = lM / lMo

    # ------------------------------------------------------------------ #
    # Active force-length characteristic  FMltilde
    # ------------------------------------------------------------------ #
    b11, b21, b31, b41 = Faparam[0], Faparam[1], Faparam[2], Faparam[3]
    b12, b22, b32, b42 = Faparam[4], Faparam[5], Faparam[6], Faparam[7]
    b13, b23, b33, b43 = 0.1, 1.0, 0.5 * _sqrt(0.5), 0.0

    FMtilde1 = b11 * _exp(-0.5 * ((lMtilde - b21) / (b31 + b41 * lMtilde))**2)
    FMtilde2 = b12 * _exp(-0.5 * ((lMtilde - b22) / (b32 + b42 * lMtilde))**2)
    FMtilde3 = b13 * _exp(-0.5 * ((lMtilde - b23) / (b33 + b43 * lMtilde + 1e-12))**2)

    FMltilde = FMtilde1 + FMtilde2 + FMtilde3
    Fiso = FMltilde

    # ------------------------------------------------------------------ #
    # Active force-velocity characteristic  FMvtilde
    # ------------------------------------------------------------------ #
    vT      = lTs * dfse / (7.0 * _exp(35.0 * (lTtilde - 0.995)))
    cos_alpha = (lMT - lTs * lTtilde) / lM
    vM      = (vMT - vT) * cos_alpha
    vMtilde = vM / vMmax

    e1, e2, e3, e4 = Fvparam[0], Fvparam[1], Fvparam[2], Fvparam[3]
    inner = e2 * vMtilde + e3
    FMvtilde = e1 * _log(inner + _sqrt(inner**2 + 1.0)) + e4

    # ------------------------------------------------------------------ #
    # Active (contractile element) force
    # ------------------------------------------------------------------ #
    d = 0.01  # damping coefficient
    Fcetilde = a * FMltilde * FMvtilde + d * vMtilde
    Fce = FMo * Fcetilde

    # ------------------------------------------------------------------ #
    # Passive force-length characteristic
    # ------------------------------------------------------------------ #
    e0  = 0.6
    kpe = 4.0
    t5  = _exp(kpe * (lMtilde - 1.0) / e0)
    Fpetilde = ((t5 - 1.0) - Fpparam[0]) / Fpparam[1]

    # ------------------------------------------------------------------ #
    # Tendon force
    # ------------------------------------------------------------------ #
    FT = fse * FMo

    # ------------------------------------------------------------------ #
    # Hill equilibrium error  (dimensionless form)
    # ------------------------------------------------------------------ #
    err = (Fcetilde + Fpetilde) * cos_alpha - fse

    return err, FT, Fce, Fiso, vMmax, Fpetilde, lMtilde


# ---------------------------------------------------------------------------
# Thin wrappers so the code is compatible with both NumPy and CasADi symbols
# ---------------------------------------------------------------------------
def _exp(x):
    try:
        import casadi as ca
        if isinstance(x, (ca.SX, ca.MX, ca.DM)):
            return ca.exp(x)
    except ImportError:
        pass
    return np.exp(x)


def _log(x):
    try:
        import casadi as ca
        if isinstance(x, (ca.SX, ca.MX, ca.DM)):
            return ca.log(x)
    except ImportError:
        pass
    return np.log(x)


def _sqrt(x):
    try:
        import casadi as ca
        if isinstance(x, (ca.SX, ca.MX, ca.DM)):
            return ca.sqrt(x)
    except ImportError:
        pass
    return np.sqrt(x)


def _sin(x):
    try:
        import casadi as ca
        if isinstance(x, (ca.SX, ca.MX, ca.DM)):
            return ca.sin(x)
    except ImportError:
        pass
    return np.sin(x)
