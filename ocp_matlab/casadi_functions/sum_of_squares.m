function J_func = sum_of_squares(var_name, D)
% This function computes the symbolic representation of the sum of square of
% an array.

% Iputs:
%   - var (str): the name that is given to the variable for CasADi to
%   track it in the expression graph.

%   - D (int): the dimensionality of the array

import casadi.*
var = SX.sym(var_name, D);

J = sum(var.^2);

J_func = Function(['J_', var_name], {var},{J});

end