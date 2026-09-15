function guess = torqueGuess(num_actuators, time_mesh, time_col)
% --------------------------------------------------------------------------
% torqueGuess
%   This function set the inital guess for the idealised torque actuators
%   and their control variables.

% INPUTs:
%   - num_actuators (int): number of actuators
%   - time_mesh (double): vector containing the mesh end-points time values.
%   - time_col (double): vector containing the collocation points time values.

 

T_mesh = size(time_mesh, 2);
T_col = size(time_col, 1);

% Torque actuators
guess.a_a = 0.1 * ones(T_mesh, num_actuators);
guess.e_a = 0.1 * ones(T_mesh, num_actuators);
guess.a_a_col = 0.1 * ones(T_col, num_actuators);

% look for invalid numerical values
validate_guess(guess);


end


