function [x_all, acc_all] = apply_constraints(x, acc, q_all_list, independent_indices)

% This function applies the registered kinematic constraints in the
% config.json file. Furthermore, it outputs the new complete set of q,
% q_dot and the acceleration in the order expected by the OpenSim model.

% Inputs:

%  - x (MX: CasADi symbolic matrix): interwined of q and q_dot
%  - acc (MX): accelerations associated to q
%  - config_struct (struct): structure with the information on the
%    kinematic constraints
%  - q_list (cell array): list of the names of the independent coordinates coming from
%    the IK results.
%  - q_all_list (cell array): list of the all the names of the coordinates coming from
%    the OpenSim model (independent + dependent).

import casadi.*

% add kinematic coupling functions
path_to_config = fullfile("kinematic_coupling", "config.json");

% read config file
json_str = fileread(path_to_config);
config_struct = jsondecode(json_str);

q = x(1:2:end);
q_dot = x(2:2:end);

num_q = size(q_all_list, 2);

% new symbolic matrix comprising all q's: dependent and independent
q_all = MX.sym('q_all', num_q);
q_dot_all = MX.sym('q_dot_all', num_q);
acc_all = MX.sym('acc_all', num_q);

% read constraint config file
constraint_names = fieldnames(config_struct);
num_constraints = size(constraint_names, 1);

% assign independent coordinate values to new full matrix
q_all(independent_indices) = q;
q_dot_all(independent_indices) = q_dot;
acc_all(independent_indices) = acc;

% loop through constraints
for i = 1:num_constraints
    constraint_name = constraint_names{i};
    
    dependent_coord_name = config_struct.(constraint_name).("dependent_coordinate");
    independent_coord_name = config_struct.(constraint_name).("independent_coordinate");
    fprintf('q ind name: %s\n', independent_coord_name);
    fprintf('q dep name: %s\n', dependent_coord_name);
    
    coupling_type = config_struct.(constraint_name).("coupling");

    independent_coord_idx_all = find(strcmp(q_all_list, independent_coord_name));
    dependent_coord_idx_all = find(strcmp(q_all_list, dependent_coord_name));
    fprintf('q ind idx: %i\n', independent_coord_idx_all);
    fprintf('q dep idx: %i\n', dependent_coord_idx_all);
    
    independent_coord_pos = q_all(independent_coord_idx_all);
    independent_coord_vel = q_dot_all(independent_coord_idx_all);
    independent_coord_acc = acc_all(independent_coord_idx_all);

    switch coupling_type
        case 'linear'
            a = config_struct.(constraint_name).("coefficients")(1);
            b = config_struct.(constraint_name).("coefficients")(2);

            value = linear(independent_coord_pos, independent_coord_vel, ...
                independent_coord_acc, a, b);

        case 'cubic'
            a = config_struct.(constraint_name).("coefficients")(1);
            b = config_struct.(constraint_name).("coefficients")(2);
            c = config_struct.(constraint_name).("coefficients")(3);
            d = config_struct.(constraint_name).("coefficients")(4);

            value = cubic(independent_coord_pos, independent_coord_vel, ...
                independent_coord_acc, a, b, c, d);

    end
    
    % write value into q_all
    q_all(dependent_coord_idx_all, :) = value.pos;
    q_dot_all(dependent_coord_idx_all, :) = value.vel;
    acc_all(dependent_coord_idx_all, :) = value.acc;
    fprintf('\n');

end


% now we need to re-order q and q_dot in the way OpenSim expects
x_all = reshape([q_all, q_dot_all]', num_q * 2, 1);

end