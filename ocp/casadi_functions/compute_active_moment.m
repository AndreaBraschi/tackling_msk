function M_func = compute_active_moment(D)

r = SX.sym('r', D);
force = SX.sym('force', D);

M = sum(r * force);    

M_func = Function('M', {r, force}, {M});


end