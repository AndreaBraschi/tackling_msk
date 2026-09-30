import casadi as ca
import numpy as np


def solve_NLP(opti, optionssol):
    # ------------------------------------------------------------------ #
    #  Extract guess
    # ------------------------------------------------------------------ #
    X = opti.x
    guess = opti.debug.value(X, opti.initial())

    # ------------------------------------------------------------------ #
    #  Sparsity pattern of constraint Jacobian
    # ------------------------------------------------------------------ #
    jac_sp = ca.jacobian(opti.g, X).sparsity()
    sp = ca.DM(jac_sp, 1)  # binary sparse matrix (nG x nX)

    sp_np = np.array(ca.DM.full(sp))  # dense numpy for indexing

    # Find constraints that depend on exactly 1 variable
    is_single = sp_np.sum(axis=1) == 1  # shape (nG,)

    # Find constraints with at most linear dependency (order 2 = quadratic -> False means linear/constant)
    is_linear = ~np.array(
        ca.which_depends(opti.g, X, 2, True)
    ).flatten().astype(bool)  # shape (nG,)

    # Simple bounds: single variable AND linear
    is_simple = is_single & is_linear

    # Which variable index does each simple constraint refer to?
    simple_rows = np.where(is_simple)[0]
    col = np.array([np.where(sp_np[r, :])[0][0] for r in simple_rows])  # variable indices

    # ------------------------------------------------------------------ #
    #  Read constraint bounds
    # ------------------------------------------------------------------ #
    lbg = np.array(opti.lbg).flatten()
    ubg = np.array(opti.ubg).flatten()

    # Detect f2(p)*x + f1(p): handles scaled variables correctly
    gf = ca.Function(
        'gf',
        [opti.x, opti.p],
        [
            opti.g[simple_rows],
            ca.jtimes(opti.g[simple_rows], opti.x, ca.DM.ones(opti.nx, 1))
        ]
    )
    f1, f2 = gf(0, opti.p)
    f1 = np.array(ca.evalf(f1)).flatten()
    f2 = np.array(ca.evalf(f2)).flatten()

    lb = (lbg[simple_rows] - f1) / np.abs(f2)
    ub = (ubg[simple_rows] - f1) / np.abs(f2)

    # ------------------------------------------------------------------ #
    #  Build lbx / ubx
    # ------------------------------------------------------------------ #
    lbx = np.full(opti.nx, -np.inf)
    ubx = np.full(opti.nx, np.inf)

    for i, j in enumerate(col):
        lbx[j] = max(lbx[j], lb[i])
        ubx[j] = min(ubx[j], ub[i])

    # ------------------------------------------------------------------ #
    #  Build reduced constraint vector (non-simple constraints only)
    # ------------------------------------------------------------------ #
    non_simple_rows = np.where(~is_simple)[0]
    new_g = opti.g[non_simple_rows]
    llb = lbg[non_simple_rows]
    uub = ubg[non_simple_rows]

    # ------------------------------------------------------------------ #
    #  Evaluate cost at initial guess
    # ------------------------------------------------------------------ #
    f_check = ca.Function('f_check', [opti.x], [opti.f])
    cost_at_guess = float(ca.evalf(f_check(guess)))
    print(f"Cost at initial guess (IPOPT x0): {cost_at_guess:.6f}")
    print(f"guess min: {np.min(guess):.6f}, max: {np.max(guess):.6f}")
    print(f"guess has NaN: {int(np.any(np.isnan(guess)))}, "
          f"has Inf: {int(np.any(np.isinf(guess)))}")

    # ------------------------------------------------------------------ #
    #  Build NLP and solve
    # ------------------------------------------------------------------ #
    prob = {'f': opti.f, 'x': opti.x, 'g': new_g}
    solver = ca.nlpsol('solver', 'ipopt', prob, optionssol)

    sol = solver(x0=guess, lbx=lbx, ubx=ubx, lbg=llb, ubg=uub)
    w_opt = np.array(ca.DM.full(sol['x'])).flatten()
    stats = solver.stats()
    g_opt = np.array(ca.DM.full(sol['g'])).flatten()
    lambda_g = np.array(ca.DM.full(sol['lam_g'])).flatten()
    lambda_x = np.array(ca.DM.full(sol['lam_x'])).flatten()

    return w_opt, stats, g_opt, lambda_g, lambda_x, solver