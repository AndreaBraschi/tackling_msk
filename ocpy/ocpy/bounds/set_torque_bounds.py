import numpy as np
from dataclasses import dataclass

def set_torque_bounds(num_actuators: int):
    # % This function set the bounds for the idealised torque actuators and their controls.
    #The bounds are scaled in such a way that they represents the +- 100%

    # Inputs:
    #  - num_act (int): number of idealised actuators

    @dataclass
    class bounds:
        a_a_lower: np.ndarray
        a_a_upper: np.ndarray
        e_a_lower: np.ndarray
        e_a_upper: np.ndarray

    @dataclass
    class scaling:
        a_a: int
        e_a: int


    bounds.a_a_lower = -np.ones(num_actuators)
    bounds.a_a_upper =  np.ones(num_actuators)

    # fixed scaling factor
    scaling.a_a = 1

    # controls to torque actuators (excitation)
    bounds.e_a_lower = -np.ones(num_actuators)
    bounds.e_a_upper =  np.ones(num_actuators)

    # fixed scaling factor
    scaling.e_a = 1

    return bounds, scaling

