import numpy as np
from typing import Dict

def generate_muscle_guess(num_muscles: int, time_mesh: np.ndarray, time_col: np.ndarray, scaling: Dict) -> Dict[str, np.ndarray]:
    """
    This function set the initial guess for the muscle state (activation and MTU force) and controls
    (activation and MTU force time derivatives).

    :param num_muscles: number of muscles of the OpenSim model
    :param time_mesh: vector containing the mesh end-points time values.
    :param time_col: vector containing the collocation points time values.
    :return:
    """

    guess = {}

    # number of segments
    T_mesh = len(time_mesh)
    # number of collocation points
    T_col = len(time_col)

    # muscle activation
    a_m = 0.1 * np.ones((num_muscles, T_mesh))
    a_m_col = 0.1 * np.ones((num_muscles, T_col))

    # MTU force
    mtu_f = 0.1 * np.ones((num_muscles, T_mesh))
    mtu_f_col = 0.1 * np.ones((num_muscles, T_col))


    # muscle activation time derivative
    va = 0.01 * np.ones((num_muscles, T_mesh))

    # MTU force time derivative
    mtu_df = 0.01 * np.ones((num_muscles, T_mesh))
    mtu_df_col = 0.01 * np.ones((num_muscles, T_col))

    guess['a_m'] = a_m
    guess['a_m_col'] = a_m_col
    guess['mtu_f'] = mtu_f / scaling['mtu_f']
    guess['mtu_f_col'] = mtu_f_col / scaling['mtu_f']
    guess['va'] = va / scaling['va']
    guess['mtu_df'] = mtu_df / scaling['mtu_df']
    guess['mtu_df_col'] = mtu_df_col / scaling['mtu_df']


    return guess