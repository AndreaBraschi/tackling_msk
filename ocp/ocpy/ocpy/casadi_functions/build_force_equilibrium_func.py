"""
build_force_equilibrium_func.py

Builds a CasADi Function that encodes the Hill-type muscle force
equilibrium (De Groote et al. 2016).
"""

import json
from pathlib import Path

import casadi as ca
import numpy as np
import sys, os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from muscle_model.force_equilibrium import force_equilibrium


def build_force_equilibrium_func(
    path_muscle_model: str | Path,
    num_muscles: int,
    MT_parameters: np.ndarray,
) -> ca.Function:
    """
    Build a CasADi Function for the Hill-equilibrium.

    Parameters
    ----------
    path_muscle_model : str or Path
        Directory containing the muscle model config.json.
    num_muscles : int
        Number of muscles.
    MT_parameters : np.ndarray, shape (5, num_muscles)
        Muscle-tendon parameters (see get_mt_parameters.py).

    Returns
    -------
    ca.Function
        'force_equilibrium'  →  (a, FTtilde, dFTtilde, lMT, vMT)
                                 ↓  ↓  ↓   ↓   ↓
                               (Hilldiff, FT, Fce, Fiso, vMmax)
    """
    path_muscle_model = Path(path_muscle_model)

    # Read force-length-velocity curve parameters from config
    with open(path_muscle_model / "config.json") as fh:
        config = json.load(fh)

    Fvparam = np.array(config["force_velocity_length"]["Fvparam"])
    Fpparam = np.array(config["force_velocity_length"]["Fpparam"])
    Faparam = np.array(config["force_velocity_length"]["Faparam"])

    # Symbolic inputs
    FTtilde  = ca.SX.sym("FTtilde",  num_muscles)   # normalised tendon forces
    a        = ca.SX.sym("a",        num_muscles)   # muscle activations
    dFTtilde = ca.SX.sym("dFTtilde", num_muscles)   # d(FTtilde)/dt
    lMT      = ca.SX.sym("lMT",      num_muscles)   # MT lengths
    vMT      = ca.SX.sym("vMT",      num_muscles)   # MT velocities

    # Symbolic outputs  (pre-allocate)
    Hilldiff = ca.SX(num_muscles, 1)
    FT       = ca.SX(num_muscles, 1)
    Fce      = ca.SX(num_muscles, 1)
    Fiso     = ca.SX(num_muscles, 1)
    vMmax    = ca.SX(num_muscles, 1)

    for m in range(num_muscles):
        # Slice single-muscle parameters: shape (5, 1)
        params_m = MT_parameters[:, m : m + 1]

        (
            Hilldiff[m],
            FT[m],
            Fce[m],
            Fiso[m],
            vMmax[m],
            _,
            _,
        ) = force_equilibrium(
            a[m], FTtilde[m], dFTtilde[m],
            lMT[m], vMT[m],
            params_m, Fvparam, Fpparam, Faparam,
        )

    return ca.Function(
        "force_equilibrium",
        [a, FTtilde, dFTtilde, lMT, vMT],
        [Hilldiff, FT, Fce, Fiso, vMmax],
    )
