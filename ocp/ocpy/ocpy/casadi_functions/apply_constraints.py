"""
apply_constraints.py

Applies the kinematic coupling constraints defined in
kinematic_coupling/config.json, then re-orders the full state vector
[q, q_dot] in the interleaved format expected by OpenSim/Simbody:
    x_all = [q_1, q_dot_1, q_2, q_dot_2, ...]
"""

import json
from pathlib import Path

import casadi as ca
import sys, os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from kinematic_coupling.linear import linear
from kinematic_coupling.cubic  import cubic


def apply_constraints(
    x,                          # ca.MX, interleaved [q; q_dot] of independent coords
    acc,                        # ca.MX, accelerations of independent coords
    q_all_list: list[str],      # names of ALL coordinates (independent + dependent)
    independent_indices,        # integer array / list of independent coord indices (0-based)
    path_to_config: str | Path = "kinematic_coupling/config.json",
):
    """
    Apply kinematic constraints and return the full interleaved state vector.

    Parameters
    ----------
    x : ca.MX, shape (2 * num_q_ind,)
        Interleaved independent [q, q_dot]: x[0::2] = q, x[1::2] = q_dot.
    acc : ca.MX, shape (num_q_ind,)
        Accelerations corresponding to the independent coordinates.
    q_all_list : list[str]
        Names of ALL coordinates in the OpenSim model (independent + dependent).
    independent_indices : array-like of int
        0-based indices into q_all_list that are *independent*.
    path_to_config : str or Path
        Path to the kinematic coupling config JSON.

    Returns
    -------
    x_all : ca.MX, shape (2 * num_q_all,)
        Full interleaved state vector for ALL coordinates.
    acc_all : ca.MX, shape (num_q_all,)
        Full acceleration vector for ALL coordinates.
    """
    path_to_config = Path(path_to_config)
    with open(path_to_config) as fh:
        config = json.load(fh)

    # Extract q and q_dot from the interleaved independent state
    q     = x[0::2]    # every even entry
    q_dot = x[1::2]    # every odd entry

    num_q_all = len(q_all_list)

    # Allocate full vectors
    q_all     = ca.MX.zeros(num_q_all, 1)
    q_dot_all = ca.MX.zeros(num_q_all, 1)
    acc_all   = ca.MX.zeros(num_q_all, 1)

    # Fill independent coordinates first
    ind = list(independent_indices)
    q_all[ind]     = q
    q_dot_all[ind] = q_dot
    acc_all[ind]   = acc

    # Apply constraints for each coupling
    for constraint_name, cfg in config.items():
        dep_name  = cfg["dependent_coordinate"]
        ind_name  = cfg["independent_coordinate"]
        coupling  = cfg["coupling"]
        coeffs    = cfg["coefficients"]

        dep_idx = q_all_list.index(dep_name)
        ind_idx = q_all_list.index(ind_name)

        pos = q_all[ind_idx]
        vel = q_dot_all[ind_idx]
        acc_in = acc_all[ind_idx]

        if coupling == "linear":
            a, b = coeffs[0], coeffs[1]
            val = linear(pos, vel, acc_in, a, b)

        elif coupling == "cubic":
            a, b, c, d = coeffs[0], coeffs[1], coeffs[2], coeffs[3]
            val = cubic(pos, vel, acc_in, a, b, c, d)

        else:
            raise ValueError(f"Unknown coupling type: '{coupling}'")

        q_all[dep_idx]     = val["pos"]
        q_dot_all[dep_idx] = val["vel"]
        acc_all[dep_idx]   = val["acc"]

    # Re-interleave into the format OpenSim/Simbody expects:
    # x_all = [q_1, q_dot_1, q_2, q_dot_2, ...]
    x_all = ca.reshape(
        ca.horzcat(q_all, q_dot_all).T,   # (2, num_q_all)
        (num_q_all * 2, 1),
    )

    return x_all, acc_all
