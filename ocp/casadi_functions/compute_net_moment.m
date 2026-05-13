function M_func = compute_net_moment(var_name, D)

r = SX.sym('r', D);
force = SX.sym('force', D);

M = sum(r * force);    

M_func = Function(['M_', var_name], {r, force}, {M});


end