"""
linear.py

Linear kinematic coupling:
    pos  = a * x + b
    vel  = a * x_vel
    acc  = a * x_acc

Works with both plain NumPy arrays and CasADi MX/SX symbolics.
"""


def linear(x, x_vel, x_acc, a: float, b: float) -> dict:
    """
    Compute position, velocity and acceleration for a linear coupling.

    Parameters
    ----------
    x     : scalar / array / CasADi symbol – independent coordinate position.
    x_vel : same type – independent coordinate velocity.
    x_acc : same type – independent coordinate acceleration.
    a, b  : float – linear coefficients  (pos = a*x + b).

    Returns
    -------
    dict with keys 'pos', 'vel', 'acc'.
    """
    y = {}
    y["pos"] = a * x + b
    y["vel"] = a * x_vel
    y["acc"] = a * x_acc
    return y
