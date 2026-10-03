import numpy as np
from scipy.interpolate import CubicSpline
from typing import Optional
from dataclasses import dataclass
# local imports
from ..utils.eval_spline import eval_spline

def set_x_bounds(Qs: np.ndarray, time_vec: np.ndarray, q_ind_indices: np.ndarray[int], num_q: int,
                 special_q_bounds):
    """
     This function set the bounds to X, which is the vector representing q and q_dot concatenated.
    :param Qs:
    :param q_ind_indices:
    :param num_q:
    :param special_q_bounds:
    :return:
    """
    @dataclass
    class q_bounds:
        X_lower: np.ndarray
        X_upper: np.ndarray
        Qsdotdot_lower: np.ndarray
        Qsdotdot_upper: np.ndarray

    @dataclass
    class scaling:
        Qs: np.ndarray
        Qsdot: np.ndarray
        Qsdotdot: np.ndarray

    # First, we approximate 1st and 2nd derivative of Qs using analytical cubic spline derivation.
    # loop through the columns of Qs, starting from index 2 (1 is time), and calculate the spline coefficients
    # from the experimental data.
    cs = CubicSpline(time_vec, Qs)
    pos, vel, acc = eval_spline(cs, time_vec, 2)


    # for the bounds, we follow the same index system of the coordinateSet object, which reflects the order of the
    # coordinates that were assigned to the .osim model. We can just loop through the coordinateSet object.
    # ----------- positions ------------ #
    Qs_lower = np.min(pos, axis=0)     # along time
    Qs_upper = np.max(pos, axis=0)     # along time

    # now let's find the coordinates that were zeroed and assign them [-1, 1] bounds
    zero_indices = np.where((Qs_lower == 0) & (Qs_upper == 0))[0]
    Qs_lower[zero_indices] = -1
    Qs_upper[zero_indices] = 1

    # scale bounds
    scaling.Qs = np.maximum(np.abs(Qs_lower), np.abs(Qs_upper))
    Qs_lower = Qs_lower / scaling.Qs
    Qs_upper = Qs_upper / scaling.Qs

    # check if the user has set some special bounds to be assigned
    if special_q_bounds:
        keys = special_q_bounds.keys()

        for key in list(keys):
            items = special_q_bounds[key]
            bounds = items["values"]
            indices = items["indices"]

            Qs_lower[indices] = bounds[0]
            Qs_upper[indices] = bounds[1]

    Qs_lower_ind = Qs_lower[q_ind_indices]
    Qs_upper_ind = Qs_upper[q_ind_indices]


    # ----------- velocities ----------- #
    Qsdot_lower = np.min(vel, axis=0)     # along time
    Qsdot_upper = np.max(vel, axis=0)

    zero_indices = np.where((Qsdot_lower == 0) & (Qsdot_upper == 0))[0]
    Qsdot_lower[zero_indices] = -1
    Qsdot_upper[zero_indices] = 1

    # scale
    scaling.Qsdot = np.maximum(np.abs(Qsdot_lower), np.abs(Qsdot_upper))
    Qsdot_lower = Qsdot_lower / scaling.Qsdot
    Qsdot_upper = Qsdot_upper / scaling.Qsdot

    Qsdot_lower_ind = Qsdot_lower[q_ind_indices]
    Qsdot_upper_ind = Qsdot_upper[q_ind_indices]


    # ----------- accelerations ----------- #
    Qsdotdot_lower = np.min(acc, axis=0)     # along time
    Qsdotdot_upper = np.max(acc, axis=0)

    zero_indices = np.where((Qsdotdot_lower == 0) & (Qsdotdot_upper == 0))[0]
    Qsdotdot_lower[zero_indices] = -1
    Qsdotdot_upper[zero_indices] = 1

    # scale
    scaling.Qsdotdot = np.maximum(np.abs(Qsdotdot_lower), np.abs(Qsdotdot_upper))
    Qsdotdot_lower = Qsdotdot_lower / scaling.Qsdotdot
    Qsdotdot_upper = Qsdotdot_upper / scaling.Qsdotdot

    q_bounds.Qsdotdot_lower_ind = Qsdotdot_lower[q_ind_indices]
    q_bounds.Qsdotdot_upper_ind = Qsdotdot_upper[q_ind_indices]

    # Now, the way Simbody/Opensim expect the q-part of the state vector isn't simply [q, q_dot],
    # but the individual dimensions of q and q_dot are rather intertwined as follows:
    # Q = [q(:, 1), q_dot(:, 1), q(:, 2), q_dot(:, 2), ...]
    # Therefore, we need to make sure that the bounds follow the same pattern,
    # as they will be assigned to the X design variables!
    X_lower = np.stack([Qs_lower_ind, Qsdot_lower_ind], axis=-1)
    X_lower = X_lower.flatten()

    X_upper = np.stack([Qs_upper_ind, Qsdot_upper_ind], axis=-1)
    X_upper = X_upper.flatten()

    q_bounds.X_lower = X_lower
    q_bounds.X_upper = X_upper


    return q_bounds, scaling
