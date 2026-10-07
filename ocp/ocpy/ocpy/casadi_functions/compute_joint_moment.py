import casadi as ca


def compute_joint_moment(D: int) -> ca.Function:

    r  = ca.SX.sym("r", D)
    force = ca.SX.sym("force", D)

    M = ca.sum1(r * force)

    return ca.Function("M", [r, force], [M])
