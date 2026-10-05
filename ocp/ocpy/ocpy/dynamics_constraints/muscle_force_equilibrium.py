import numpy as np


def force_equilibrium(a_m, fT_norm, dfT_norm, lMT, vMT, mtu_params):
    """
    Compute the Hill-equilibrium and related muscle quantities.

    """

    # ---------------------------------------------------------------------------------------------------------------- #
    #                                        static parameters to form muscle ODEs
    # ---------------------------------------------------------------------------------------------------------------- #
    # Unpack muscle-tendon parameters (broadcast over muscles)
    FMo   = mtu_params[0, :]        # peak isometric muscle force
    lMo   = mtu_params[1, :]        # optimal fiber length
    lTs   = mtu_params[2, :]        # tendon slack length
    alpha_o = mtu_params[3, :]       # pennation angle at optimal fiber length
    vM_max  = mtu_params[4, :]       # maximal muscle fiber velocity


    # Tendon force-length
    kT = 35.0
    c1 = 0.2
    c2 = 0.995
    c3 = 0.25


    # static parameters for the muscle-length-velocity relationships (active and passive)
    f_p_params = [-0.9952, 53.5982]                                                            # passive force-length
    f_v_params = [-0.3183, -8.1492, -0.3741, 0.8856]                                           # force-velocity

    f_a_params = [0.8145, 1.0550, 0.1624, 0.0633, 0.4330, 0.7168, -0.0299, 0.2004]  # active force-length
    b11, b12, b13 = f_a_params[0], f_a_params[4], 0.1
    b21, b22, b23 = f_a_params[1], f_a_params[5], 1.0
    b31, b32, b33 = f_a_params[2], f_a_params[6], 0.354
    b41, b42, b43 = f_a_params[3], f_a_params[7], 0.0

    # ------------------------------------------------------------------ #
    # Tendon length
    # ------------------------------------------------------------------ #
    # Invert tendon force-length: recover normalised tendon length
    lT_norm = _log((fT_norm + c3) / c1) / kT + c2
    # tendon length (unscaled)
    lT = lT_norm * lTs

    # ------------------------------------------------------------------ #
    # Fibre length
    # ------------------------------------------------------------------ #
    lM = _sqrt((lMo * _sin(alpha_o))**2 + (lMT - lT)**2)         # muscle fiber length
    lM_norm = lM / lMo

    # ------------------------------------------------------------------ #
    # Pennation angle
    # ------------------------------------------------------------------ #
    cos_alpha = (lMT - lT) / lM

    # ------------------------------------------------------------------ #
    # Normalised Tendon velocity
    # ------------------------------------------------------------------ #
    vT_norm = dfT_norm / (c1 * kT * _exp(kT * (lT_norm - c2)))
    vT = vT_norm * lTs

    # ------------------------------------------------------------------ #
    # Muscle fibers velocity
    # ------------------------------------------------------------------ #
    vM = (vMT - vT) * cos_alpha
    vM_norm = vM / vM_max


    # Now to compute the Hill-equilibrium, we need the following quantities: active force-length, passive force-length
    # and force-velocity.

    # ------------------------------------------------------------------ #
    # Active force-length function
    # ------------------------------------------------------------------ #
    # 3 Gaussian functions
    g_func_1 = b11 * _exp(-0.5 * ((lM_norm - b21) / (b31 + b41 * lM_norm))**2)
    g_func_2 = b12 * _exp(-0.5 * ((lM_norm - b22) / (b32 + b42 * lM_norm))**2)
    g_func_3 = b13 * _exp(-0.5 * ((lM_norm - b23) / (b33 + b43 * lM_norm + 1e-12))**2)

    fl_act = g_func_1 + g_func_2 + g_func_3


    # ------------------------------------------------------------------ #
    # Passive force-length function
    # ------------------------------------------------------------------ #
    k_pe = 4.0
    e0 = 0.6
    exp_0 = _exp((k_pe * (lM_norm - 1)) / e0) - 1
    exp_1 = _exp(k_pe) - 1
    fl_passive = (exp_0 - f_p_params[0]) / exp_1

    # ------------------------------------------------------------------ #
    # Force-Velocity function
    # ------------------------------------------------------------------ #
    d0, d1, d2, d3 = f_v_params[0], f_v_params[1], f_v_params[2], f_v_params[3]
    fv_term = d1 * vM_norm + d2
    fv = d0 * _log((fv_term + _sqrt(fv_term ** 2 + 1))) + d3


    # ------------------------------------------------------------------ #
    # Hill equilibrium error  (dimensionless form)
    # ------------------------------------------------------------------ #
    muscle_force = a_m * fl_act * fv + fl_passive
    err = muscle_force * cos_alpha - fT_norm

    return err


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
