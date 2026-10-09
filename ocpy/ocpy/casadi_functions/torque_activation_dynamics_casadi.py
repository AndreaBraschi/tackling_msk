"""
torque_activation_dynamics_casadi.py

Wraps the torque activation dynamics ODE in a CasADi Function so it can
participate in the CasADi expression graph / AD machinery.
"""
import casadi as ca
from ..dynamics_constraints.torque_activation import torque_activation_dynamics

def torque_activation_dynamics_casadi(num_actuators: int) -> ca.Function:
    """
    Build a CasADi Function for the torque actuator activation dynamics.

    Parameters
    ----------
    num_actuators : Number of torque actuators of the model.

    Returns
    -------
    ca.Function
        Maps (excitation, activation) → da/dt.
        Name: 'torque_activation_dynamics_casadi'.
    """
    excitation = ca.SX.sym("excitation", num_actuators)
    activation = ca.SX.sym("activation", num_actuators)

    dadt = torque_activation_dynamics(excitation, activation)

    return ca.Function(
        "torque_activation_dynamics_casadi",
        [excitation, activation],
        [dadt],
    )
