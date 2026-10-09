import numpy as np
import casadi as ca

def radau(d: int, method: str = "radau"):
    """Return (tau_root, C, D, B) for a d-point collocation scheme."""
    tau_root = np.concatenate([[0], ca.collocation_points(d, method)])
    C = np.zeros((d + 1, d + 1))
    D = np.zeros(d + 1)
    B = np.zeros(d + 1)

    for j in range(d + 1):
        p = np.poly1d([1])
        for r in range(d + 1):
            if r != j:
                p *= np.poly1d([1, -tau_root[r]]) / (tau_root[j] - tau_root[r])
        D[j] = np.polyval(p, 1.0)
        dp = np.polyder(p)
        for r in range(d + 1):
            C[j, r] = np.polyval(dp, tau_root[r])
        pint = np.polyint(p)
        B[j] = np.polyval(pint, 1.0)

    return tau_root, C, D, B