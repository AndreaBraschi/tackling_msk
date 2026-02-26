from math import sqrt


def compute_fp(x: float) -> float:

    """
    :param x: penetration between the 2 geometries. This is a scalar value. It's a 3D vector projected onto the normal
              unit vector between the 2 geometries.
    :return:
    """

    f_p: float = x ** 1.5
    return f_p

def compute_fv(x_dot: float, c: float) -> float:
    """
    :param x_dot: rate of change of the penetration. This is also a scalar value as described above.
    :param c: damping.
    :return:
    """
    f_v: float = 1 + (1.5 * c * x_dot)
    return f_v


def hunt_crossley(f_p: float, f_v: float, sphere_r: float, k: float) -> float:
    """
    This function computes the Hunt-Crossley contact force.

    :param sphere_r:  radius of sphere.
    :param k: stiffness.
    :return: force: scalar value that is to be projected on the same normal vector used to compute the scalar penetration
                    and its rate of change over time.
    """

    f_n: float = (4/3) * (k ** 1.5) * sqrt(sphere_r) * f_p * f_v

    return f_n