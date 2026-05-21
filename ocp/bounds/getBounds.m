function [bounds, scaling] = getBounds(Qs, GRFs, num_q, num_muscles, num_act, grf_indices, grm_indices)

% This function assign the bounds to: model coordinates (Qs), Ground
% Reaction Forces (GRF), muscle activation, rate of change of muscle
% activation over time, MT force and rate of change of MT force over time. 
% Then, the bounds are scaled in such a way that they represents the +- 100%

%
% Inputs:
%   - Qs (struct): 
%   - GRF (struct):
%   - num_q (int):
%   - num_muscles (int):
%   - num_act (int):
%   - grf_indices (int):

import org.opensim.modeling.*

time = Qs.allfilt(:, 1);
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
end

% ----------------------------- Qs bounds ----------------------------- %
% for the bounds, we follow the same index system of the coordinateSet
% object, which reflects the order of the coordinates that were assigned to
% the .osim model. We can just loop through the coordinateSet object.

% --- positions --- %
for i = 1:num_q
    Qs.upper(i) = max(y(i).pos);
    Qs.lower(i) = min(y(i).pos);
end
% scale
scaling.Qs = max(abs(Qs.lower), abs(Qs.upper));
Qs.lower = (Qs.lower)./scaling.Qs;
Qs.upper = (Qs.upper)./scaling.Qs;


% --- velocities --- %
for i = 1:num_q
    Qsdot.upper(i) = max(y(i).vel); 
    Qsdot.lower(i) = min(y(i).vel);
end
% scale
scaling.Qsdot = max(abs(Qsdot.lower), abs(Qsdot.upper));
Qsdot.lower = (Qsdot.lower)./scaling.Qsdot;
Qsdot.upper = (Qsdot.upper)./scaling.Qsdot;

% --- accelerations --- %
for i = 1:num_q
    bounds.Qsdotdot.upper(i) = max(y(i).acc);
    bounds.Qsdotdot.lower(i) = min(y(i).acc);
end
% scale
scaling.Qsdotdot = max(abs(bounds.Qsdotdot.lower),abs(bounds.Qsdotdot.upper));
bounds.Qsdotdot.lower = (bounds.Qsdotdot.lower)./scaling.Qsdotdot;
bounds.Qsdotdot.upper = (bounds.Qsdotdot.upper)./scaling.Qsdotdot;


% Now, the way Simbody/Opensim expect the q-part of the state vector isn't
% simply [q, q_dot], but the individual dimensions of q and q_dot are
% rather interwinded as follos:
% Q = [q(:, 1), q_dot(:, 1), q(:, 2), q_dot(:, 2), ...]
% Therefore, we need to make sure that the bounds follow the same pattern,
% as they will be assigned to the X design variables!
X_lower = reshape([Qs.lower, Qsdot.lower], 2, num_q);
X_lower = reshape(permute(X_lower, [2, 1]), 1, dims);

X_upper = reshape([Qs.upper, Qsdot.upper], 2, num_q);
X_upper = reshape(permute(X_upper, [2, 1]), 1, dims);

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


% -------------------------- Torque Actuators -------------------------- %
% Torque actuator activations
bounds.a_a.lower = -ones(1, num_act);
bounds.a_a.upper = ones(1, num_act);

% Torque actuator excitation
bounds.e_a.lower = -ones(1, num_act);
bounds.e_a.upper = ones(1, num_act);


% ------------------------- Contact Model bounds ------------------------- %

% ----> why there was no bounds for stiffness and damping? 
% ----> were they fixed?
% ----> what should I do? give full flexibility of constraining to some
% mechanically meaningful values?

end