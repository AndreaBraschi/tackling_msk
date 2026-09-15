function [bounds, scaling] = set_xBounds(Qs, independent_coord_idx, num_q, special_q_bounds)

% This function set the bounds to X, which is the vector representing q and
% q_dot concatenated.

% Inputs:
%   - Qs (struct): 
%   - independent_coord_idx (double): 
%   - num_q (int): number of independent q's
%   - special_q_bounds (struct): section of the config file where the user
%   may have set specific bounds for certain variables.
   

import org.opensim.modeling.*

time = Qs.time;
y = struct('pos', {}, 'vel', {}, 'acc', {});

% Approximate 1st and 2nd derivative of Qs using analytical cubic spline
% derivation.
% loop through the columns of Qs, starting from index 2 (1 is time), and
% calculate the spline coefficients, given the experimental data
for i = 1:num_q
    
    cs = spline(time, Qs.allfilt(:, i + 1));
    y(i) = eval_spline(cs, time, 2);


end

% ----------------------------- Qs bounds ----------------------------- %
% for the bounds, we follow the same index system of the coordinateSet
% object, which reflects the order of the coordinates that were assigned to
% the .osim model. We can just loop through the coordinateSet object.

% --- positions --- %
for i = 1:num_q
    if max(y(i).pos) == 0 || min(y(i).pos) == 0
        fprintf('WARNING: %d coordinate(s) has been detected:\n', i);
        Qs.upper(i) = 1;
        Qs.lower(i) = -1;

    else
        Qs.upper(i) = max(y(i).pos); 
        Qs.lower(i) = min(y(i).pos);
    end

    if Qs.lower(i) == Qs.upper(i)
        fprintf('WARNING: %d coordinate(s) have identical upper and lower bounds:\n', i);
    end
    

end
% scale
scaling.Qs = max(abs(Qs.lower), abs(Qs.upper));
Qs.lower = (Qs.lower)./scaling.Qs;
Qs.upper = (Qs.upper)./scaling.Qs;

% check if the user has set some special bounds to be assigned
if ~isempty(special_q_bounds)
    keys = fieldnames(special_q_bounds);
    for i = 1:length(keys)

        key   = keys{i};
        items = special_q_bounds.(key);
        special_bounds = items.values;   % values of bounds: [lower, upper]
        indices = items.indices;         % indices to where to apply the bounds

        Qs.lower(:, indices) = special_bounds(1);
        Qs.upper(:, indices) = special_bounds(2);


    end
 
end

Qs.lower_ind = Qs.lower(:, independent_coord_idx);
Qs.upper_ind = Qs.upper(:, independent_coord_idx);


% --- velocities --- %
for i = 1:num_q    
    if max(y(i).vel) == 0 || min(y(i).vel) == 0
        Qsdot.upper(i) = 1;
        Qsdot.lower(i) = -1;

    else
        Qsdot.upper(i) = max(y(i).vel); 
        Qsdot.lower(i) = min(y(i).vel);
    end

    if Qsdot.lower(i) == Qsdot.upper(i)
        fprintf('WARNING: %d coordinate(s) have identical upper and lower bounds:\n', i);
    end
end

% scale
scaling.Qsdot = max(abs(Qsdot.lower), abs(Qsdot.upper));
Qsdot.lower = (Qsdot.lower)./scaling.Qsdot;
Qsdot.upper = (Qsdot.upper)./scaling.Qsdot;


Qsdot.lower_ind = Qsdot.lower(:, independent_coord_idx');
Qsdot.upper_ind = Qsdot.upper(:, independent_coord_idx');


% --- accelerations --- %
for i = 1:num_q
    if max(y(i).acc) == 0 || min(y(i).acc) == 0
        bounds.Qsdotdot.upper(i) = 1;
        bounds.Qsdotdot.lower(i) = -1;

    else
        bounds.Qsdotdot.upper(i) = max(y(i).acc);
        bounds.Qsdotdot.lower(i) = min(y(i).acc);
    end
end

% scale
scaling.Qsdotdot = max(abs(bounds.Qsdotdot.lower), abs(bounds.Qsdotdot.upper));
bounds.Qsdotdot.lower = (bounds.Qsdotdot.lower)./scaling.Qsdotdot;
bounds.Qsdotdot.upper = (bounds.Qsdotdot.upper)./scaling.Qsdotdot;

bounds.Qsdotdot.lower_ind = bounds.Qsdotdot.lower(:, independent_coord_idx');
bounds.Qsdotdot.upper_ind = bounds.Qsdotdot.upper(:, independent_coord_idx');


% Now, the way Simbody/Opensim expect the q-part of the state vector isn't
% simply [q, q_dot], but the individual dimensions of q and q_dot are
% rather interwinded as follos:
% Q = [q(:, 1), q_dot(:, 1), q(:, 2), q_dot(:, 2), ...]
% Therefore, we need to make sure that the bounds follow the same pattern,
% as they will be assigned to the X design variables!
dims_ind = size(independent_coord_idx, 2) * 2;
X_lower = cat(3, Qs.lower_ind, Qsdot.lower_ind);
X_lower = reshape(permute(X_lower, [1, 3, 2]), 1, dims_ind);

X_upper = cat(3, Qs.upper_ind, Qsdot.upper_ind);
X_upper = reshape(permute(X_upper, [1, 3, 2]), 1, dims_ind);


bounds.X.lower = X_lower;
bounds.X.upper = X_upper;

% -------------------------- Validate struct -------------------------- %
validate_bounds(bounds);

end