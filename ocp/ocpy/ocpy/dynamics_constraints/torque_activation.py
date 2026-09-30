"""
torque_activation.py

First-order activation dynamics for torque actuators:
    da/dt = (e - a) / tau

Works with NumPy arrays and CasADi symbolics alike.
"""


def torque_activation_dynamics(e, a):
    """
    Compute the time derivative of torque actuator activation.

    Parameters
    ----------
    e : array-like or CasADi symbol
        Excitation signal(s).
    a : array-like or CasADi symbol
        Current activation(s).

    Returns
    -------
    da_dt : same type as inputs
        Time derivative of activation.
    """
    tau = 0.035
    return (e - a) / tau
