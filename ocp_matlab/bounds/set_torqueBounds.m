function [bounds, scaling] = set_torqueBounds(num_act)

% This function set the bounds for the idealised torque actuators and their
% controls. 
% The bounds are scaled in such a way that they represents the +- 100%

%
% Inputs:
%   - num_act (int): number of idealised actuators

% Torque actuator activations
bounds.a_a.lower = -ones(1, num_act);
bounds.a_a.upper = ones(1, num_act);

% fixed scaling factor
scaling.a_a = 1;

% Torque actuator excitation
bounds.e_a.lower = -ones(1, num_act);
bounds.e_a.upper = ones(1, num_act);

% fixed scaling factor
scaling.e_a = 1;


% -------------------------- Validate struct -------------------------- %
validate_bounds(bounds);

end