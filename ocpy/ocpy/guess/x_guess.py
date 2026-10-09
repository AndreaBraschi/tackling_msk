import numpy as np
from scipy.interpolate import CubicSpline
from typing import Dict

# internal imports
from ..utils.eval_spline_col import *

def generate_x_guess(Qs: np.ndarray, q_ind_indices: np.ndarray[int], time_experimental: np.ndarray,
                     time_mesh: np.ndarray, time_col: np.ndarray, num_q: int, scaling) -> Dict[str, np.ndarray]:

    """
    This function sets the initial guess for X (q's and q_dot's interleaved) from the experimental data.

    The experimental data is used to compute cubic spline coefficients, so
    that a prediction on the values of q's can be made at the various
    collocation points. Furthermore, velocity and acceleration of the
    splines are computed at the same points to form a final initial guess.

    This function was adapted to the specific case of this repo, but the
    original, generic version can be found at:
    https://github.com/KULeuvenNeuromechanics/PredSim/OCP/getGuess_DI_opti.m

    :param Qs: generalised coordinates
    :param q_ind_indices: indices of the independent generalised coordinates
    :param time_experimental: time vector from experimental data
    :param time_mesh: time vector representing the time at each first point of the collocation meshes
    :param time_col: time vector representing the time at each collocation point of the trajectory
    :param num_q: total nuber of the model's generalised coordinates
    :param scaling: scaling factor of the model's generalised coordinates
    :return:
    """

    # initialise the 'guess' object as a Dict
    guess: Dict[str, np.ndarray] = {}

    mesh_k = np.searchsorted(time_experimental, time_mesh, side='right') - 1
    mesh_dt = time_mesh - time_experimental[mesh_k]

    col_k = np.searchsorted(time_experimental, time_col, side='right') - 1
    col_dt =  time_col -  time_experimental[col_k]

    # interpolate the generalised coordinates using CubicSpline
    cs = CubicSpline(time_experimental, Qs)     # cubic spline coefficients

    try:
        pos_end, vel_end, acc_end = eval_spline_col(cs, time_mesh, mesh_k, mesh_dt[:, None])
    except:
        mesh_k = np.clip(mesh_k, 0, len(time_experimental) - 2)  # clips 120 → 119
        pos_end, vel_end, acc_end = eval_spline_col(cs, time_mesh, mesh_k, mesh_dt[:, None])


    # we might be in the situation where the the last index inside col_k is greater than the number of spline segments.
    # However, the col_k indices represent a correspondence of points and, for a given spline segment, the 2 knots
    # (first and last point of the segments) are represented by the same polynomial. Hence why we can clip the point to
    # the one before.
    try:
        pos_col, vel_col, acc_col = eval_spline_col(cs, time_col, col_k, col_dt[:, None])
    except:
        col_k = np.clip(col_k, 0, len(time_experimental) - 2)  # clips 120 → 119
        pos_col, vel_col, acc_col = eval_spline_col(cs, time_col, col_k, col_dt[:, None])


    # -------------------- scale ------------------- #
    # mesh end points
    pos_end_scaled = pos_end / scaling.Qs
    vel_end_scaled = vel_end / scaling.Qsdot
    acc_end_scaled = acc_end / scaling.Qsdotdot

    # collocation points
    pos_col_scaled = pos_col / scaling.Qs
    vel_col_scaled = vel_col / scaling.Qsdot
    acc_col_scaled = acc_col / scaling.Qsdotdot

    # Now, the way Simbody/Opensim expect the q-part of the state vector isn't simply [q, q_dot],
    # but the individual dimensions of q and q_dot are rather intertwined as follows:
    # Q = [q(:, 1), q_dot(:, 1), q(:, 2), q_dot(:, 2), ...]
    # Therefore, we need to make sure that the bounds follow the same pattern,
    # as they will be assigned to the X design variables!
    x_end = np.stack([pos_end_scaled, vel_end_scaled], axis=-1)
    x_end = x_end.reshape(*x_end.shape[:-2], -1)
    guess['x_end'] = x_end

    x_col = np.stack([pos_col_scaled, vel_col_scaled], axis=-1)
    x_col = x_col.reshape(*x_col.shape[:-2], -1)
    guess['x_col'] = x_col

    # Do the same for the independent coordinates only.
    x_end_ind = np.stack([pos_end_scaled[:, q_ind_indices], vel_end_scaled[:, q_ind_indices]], axis=-1)
    x_end_ind = x_end_ind .reshape(*x_end_ind.shape[:-2], -1)
    guess['x_end_ind'] = x_end_ind

    x_col_ind = np.stack([pos_col_scaled[:, q_ind_indices], vel_col_scaled[:, q_ind_indices]], axis=-1)
    x_col_ind = x_col_ind.reshape(*x_col_ind.shape[:-2], -1)
    guess['x_col_ind'] = x_col_ind

    # ----------- add accelerations to the dictionary ----------- #
    # all
    guess['acc_all'] = acc_col_scaled
    guess['acc_col'] = acc_col_scaled[:, q_ind_indices]

    return guess

    x = 5