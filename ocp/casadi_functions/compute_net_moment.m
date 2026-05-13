function M_func = compute_net_moment(var_name, D)

r = SX.sym(['r_', var_name], D);
force = SX.sym(['force_', var_name], D);

M = sum(r * force);    

M_func = Function(['M_', var_name], {r, force}, {M});


end