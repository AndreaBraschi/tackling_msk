function [x_all, acc_all] = apply_constraints(x, acc, q_list, q_all_list)

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
path_name = "../kinematic_coupling";
addpath(path_name);

% read config file
json_str = fileread("config.json");
config_struct = jsondecode(json_str);

q = x(1:2:end);
q_dot = x(2:2:end);

q_all = MX.sym('q_all', size(q, 1));
q_dot_all = MX.sym('q_dot_all', size(q, 1));
acc_all = MX.sym('acc_all', size(q, 1));

constraint_names = fieldnames(config_struct);
num_constraints = size(constraint_names, 1);

independent_indices = cellfun(@(q) find(strcmp(q_all_list, q)), q_list);
q_all(independent_indices) = q;
q_dot_all(independent_indices) = q_dot;
acc_all(independent_indices) = acc;

% loop through constraints
for i = 1:num_constraints
    constraint_name = constraint_names{i};
    dependent_coord_name = config_struct.(constraint_name).dependent_coordinate;
    independent_coord_name = config_struct.(constraint_name).independent_coordinate;
    coupling_type = config_struct.(constraint_name).coupling;

    independent_coord_idx = strcmp(q_list, independent_coord_name);
    dependent_coord_idx_all = strcmp(q_all_list, dependent_coord_name);
    
    independent_coord_pos = q(independent_coord_idx);
    independent_coord_vel = q_dot(independent_coord_idx);
    independent_coord_acc = acc(independent_coord_idx);

    switch coupling_type
        case 'linear'
            a = config_struct.(constraint_name).coefficients(1);
            b = config_struct.(constraint_name).coefficients(2);

            value = linear(independent_coord_pos, independent_coord_vel, independent_coord_acc, a, b);

            % write values
            q_all(dependent_coord_idx_all) = value.pos;
            q_dot_all(dependent_coord_idx_all) = value.vel;
            acc_all(dependent_coord_idx_all) = value.acc;


        case 'cubic'
            a = config_struct.(constraint_name).coefficients(1);
            b = config_struct.(constraint_name).coefficients(2);
            c = config_struct.(constraint_name).coefficients(3);
            d = config_struct.(constraint_name).coefficients(4);

            value = cubic(independent_coord_value, independent_coord_vel, independent_coord_acc, a, b, c, d);

            % write value into q_all
            q_all(dependent_coord_idx_all) = value.pos;
            q_dot_all(dependent_coord_idx_all) = value.vel;
            acc_all(dependent_coord_idx_all) = value.acc;

    end

end


dims = size(x, 1);

% now we need to re-order q and q_dot in the way OpenSim expects
x_all = reshape([q_all, q_dot_all], 2, dims);
x_all = reshape(permute(x_all, [2, 1]), 1, dims * 2);


end