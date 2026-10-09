import numpy as np
from typing import Dict

def generate_torque_guess(num_actuators: int, time_mesh: np.ndarray, time_col: np.ndarray) -> Dict[str, np.ndarray]:
    """
    This function set the inital guess for the idealised torque actuators and their control variables.

    :param num_actuators: number of actuators
    :param time_mesh: vector containing the mesh end-points time values.
    :param time_col: vector containing the collocation points time values.
    :return:
    """

    guess: Dict[str, np.ndarray] = {}

    T_mesh = len(time_mesh)
    T_col = len(time_col)


    guess['a_a'] = 0.1 * np.ones((T_mesh, num_actuators))
    guess['e_a'] = 0.1 * np.ones((T_mesh, num_actuators))
    guess['a_a_col'] = 0.1 * np.ones((T_col, num_actuators))

    return guess