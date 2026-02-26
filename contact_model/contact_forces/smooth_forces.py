from math import tanh, sqrt, pow
from .contact_forces import compute_fp, compute_fv, hunt_crossley


def smooth_hunt_crossley(x: float, x_dot: float, sphere_r: float, k: float, c: float,
                         bc: float = 50.0, cf: float = 1e-8) -> float:

    k_new: float = 0.5 * pow(k, 2 / 3)
    f_p: float = compute_fp(sqrt(x ** 2 + cf))
    f_p_smooth: float = f_p * (0.5 + 0.5 * tanh(bc * x))

    f_v: float = compute_fv(x_dot, c)
    f_v_smooth: float = f_v * (0.5 + 0.5 * tanh(bc * (x_dot + 2 / (3 * c))))

    r_new: float = sphere_r * k_new

    # now use the 2 smooth versions to compute Hunt-Crossley force
    f_hc: float = hunt_crossley(f_p_smooth, f_v_smooth, r_new, k_new ** (2 / 3))

    return f_hc


