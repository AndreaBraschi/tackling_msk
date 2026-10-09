"""
cubic.py

Cubic kinematic coupling:
    pos  = a*x^3 + b*x^2 + c*x + d
    vel  = (3*a*x^2 + 2*b*x + c) * x_vel
    acc  = (6*a*x + 2*b) * x_vel^2 + (3*a*x^2 + 2*b*x + c) * x_acc

Works with both plain NumPy arrays and CasADi MX/SX symbolics.
"""


def cubic(x, x_vel, x_acc, a: float, b: float, c: float, d: float) -> dict:
    """
    Compute position, velocity and acceleration for a cubic coupling.

    Parameters
    ----------
    x, x_vel, x_acc : scalar / array / CasADi symbol
        Independent coordinate: position, velocity, acceleration.
    a, b, c, d : float
        Cubic polynomial coefficients  (pos = a*x^3 + b*x^2 + c*x + d).

    Returns
    -------
    dict with keys 'pos', 'vel', 'acc'.
    """
    y = {}
    y["pos"] = a * x**3 + b * x**2 + c * x + d

    # first derivative w.r.t. time (chain rule)
    dy_dx = 3 * a * x**2 + 2 * b * x + c
    y["vel"] = dy_dx * x_vel

    # second derivative w.r.t. time
    d2y_dx2 = 6 * a * x + 2 * b
    y["acc"] = d2y_dx2 * x_vel**2 + dy_dx * x_acc

    return y
