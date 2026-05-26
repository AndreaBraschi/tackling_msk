function [bounds, scaling] = getBounds(Qs, independent_coord_idx, GRFs, num_q, num_muscles, num_act, grf_indices, grm_indices)

% This function assign the bounds to: model coordinates (Qs), Ground
% Reaction Forces (GRF), muscle activation, rate of change of muscle
% activation over time, MT force and rate of change of MT force over time. 
% Then, the bounds are scaled in such a way that they represents the +- 100%

%
% Inputs:
%   - Qs (struct): 
%   - independent_coord_idx (double): 
%   - GRF (struct):
%   - num_muscles (int):
%   - num_act (int):
%   - grf_indices (int):

import org.opensim.modeling.*

time = Qs.time;
T = size(time, 1);  % period

y = struct('pos', {}, 'vel', {}, 'acc', {});
dims = num_q * 2;

% Approximate 1st and 2nd derivative of Qs using analytical cubic spline
% derivation.
% loop through the columns of Qs, starting from index 2 (1 is time), and
% calculate the spline coefficients, given the experimental data
for i = 1:num_q
    
    cs = spline(time, Qs.allfilt(:, i + 1));
    y(i) = eval_spline(cs, time, 2);

    %validate_array(y(i).pos, "Q_spline")
    %validate_array(y(i).vel, "Qdots_spline")
    %validate_array(y(i).acc, "Qdotdots_spline")

end

% ----------------------------- Qs bounds ----------------------------- %
% for the bounds, we follow the same index system of the coordinateSet
% object, which reflects the order of the coordinates that were assigned to
% the .osim model. We can just loop through the coordinateSet object.

% --- positions --- %
for i = 1:num_q
    if max(y(i).pos) == 0 || min(y(i).pos) == 0
        Qs.upper(i) = 1;
        Qs.lower(i) = 1;

    else
        Qs.upper(i) = max(y(i).pos); 
        Qs.lower(i) = min(y(i).pos);
    end
end
% scale
scaling.Qs = max(abs(Qs.lower), abs(Qs.upper));
Qs.lower = (Qs.lower)./scaling.Qs;
Qs.upper = (Qs.upper)./scaling.Qs;

Qs.lower_ind = Qs.lower(:, independent_coord_idx);
Qs.upper_ind = Qs.upper(:, independent_coord_idx);


% --- velocities --- %
for i = 1:num_q    
    if max(y(i).vel) == 0 || min(y(i).vel) == 0
        Qsdot.upper(i) = 1;
        Qsdot.lower(i) = 1;

    else
        Qsdot.upper(i) = max(y(i).vel); 
        Qsdot.lower(i) = min(y(i).vel);
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
        bounds.Qsdotdot.lower(i) = 1;

    else
        bounds.Qsdotdot.upper(i) = max(y(i).acc);
        bounds.Qsdotdot.lower(i) = min(y(i).acc);
    end
end

% scale
scaling.Qsdotdot = max(abs(bounds.Qsdotdot.lower), abs(bounds.Qsdotdot.upper));
bounds.Qsdotdot.lower = (bounds.Qsdotdot.lower)./scaling.Qsdotdot;
bounds.Qsdotdot.upper = (bounds.Qsdotdot.upper)./scaling.Qsdotdot;


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

% ----------------------------- GRFs bounds ----------------------------- %
lower_grf = min(GRFs.data(:, grf_indices));
upper_grf = max(GRFs.data(:, grf_indices));

% extend bounds to give some flexibility
force_range = abs(upper_grf - lower_grf);
bounds.GRF.lower = lower_grf - force_range;
bounds.GRF.upper = upper_grf + force_range;

% scale
scaling.GRF = max(abs(bounds.GRF.lower), abs(bounds.GRF.upper));
bounds.GRF.lower = (bounds.GRF.lower)./scaling.GRF;
bounds.GRF.upper = (bounds.GRF.upper)./scaling.GRF;


% ----------------------------- GRMs bounds ----------------------------- %
lower_grm = min(GRFs.data(:, grm_indices));
upper_grm = max(GRFs.data(:, grm_indices));

% extend bounds to give some flexibility
moment_range = abs(upper_grm - lower_grm);
bounds.GRM.lower = lower_grm - moment_range;
bounds.GRM.upper = upper_grm + moment_range;

% scale
scaling.GRM = max(abs(bounds.GRM.lower), abs(bounds.GRM.upper));
bounds.GRM.lower = (bounds.GRM.lower)./scaling.GRM;
bounds.GRM.upper = (bounds.GRM.upper)./scaling.GRM;


% ---------------------------- Muscle bounds ---------------------------- %
% ---- activation ---- %
bounds.a.lower = zeros(1, num_muscles);
bounds.a.upper = ones(1, num_muscles);

% scaling
scaling.a = 1;

% ---- d_activation / dt (time derivative) ---- %
tact = 0.015;
tdeact = 0.06;
bounds.vA.lower = (-1/100 * ones(1, num_muscles))./(ones(1, num_muscles) * tdeact);
bounds.vA.upper = (1/100 * ones(1, num_muscles))./(ones(1, num_muscles) * tact);

% fixed scaling factor
scaling.vA = 100;

% ---- MT forces ---- %
bounds.FTtilde.lower = zeros(1, num_muscles);
bounds.FTtilde.upper = 5 * ones(1, num_muscles);

% scale ---> why here yes and vA no?
scaling.FTtilde  = max(abs(bounds.FTtilde.lower),abs(bounds.FTtilde.upper)); 
bounds.FTtilde.lower = (bounds.FTtilde.lower)./scaling.FTtilde;
bounds.FTtilde.upper = (bounds.FTtilde.upper)./scaling.FTtilde;

% ---- d_MTf / dt ---- %
bounds.dFTtilde.lower = -1*ones(1, num_muscles);
bounds.dFTtilde.upper = 1*ones(1, num_muscles);

% fixed scaling factor
scaling.dFTtilde = 100;

% -------------------------- Torque Actuators -------------------------- %
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

% look for invalid numerical values
validate_bounds(bounds);

% ------------------------- Contact Model bounds ------------------------- %

% ----> why there was no bounds for stiffness and damping? 
% ----> were they fixed?
% ----> what should I do? give full flexibility of constraining to some
% mechanically meaningful values?

end


function validate_array(val, name)
    nan_idx = find(isnan(val(:)));
    inf_idx = find(isinf(val(:)));
    
    if ~isempty(nan_idx)
        error('[Array Validation] "%s" contains NaN at indices: %s', ...
            name, mat2str(nan_idx'));
    end
    
    if ~isempty(inf_idx)
        error('[Array Validation] "%s" contains Inf at indices: %s', ...
            name, mat2str(inf_idx'));
    end
end


function validate_bounds(bounds)
    fields = fieldnames(bounds);
    for i = 1:numel(fields)
        field = fields{i};
        
        subfields = {'lower', 'upper'};
        for j = 1:numel(subfields)
            subfield = subfields{j};
            val = bounds.(field).(subfield);
            
            nan_idx = find(isnan(val(:)));
            inf_idx = find(isinf(val(:)));
            
            if ~isempty(nan_idx)
                error('[Bounds Validation] "%s.%s" contains NaN at indices: %s', ...
                    field, subfield, mat2str(nan_idx'));
            end
            
            if ~isempty(inf_idx)
                error('[Bounds Validation] "%s.%s" contains Inf at indices: %s', ...
                    field, subfield, mat2str(inf_idx'));
            end
        end
    end
    fprintf('[Bounds Validation] All fields passed.\n');
end