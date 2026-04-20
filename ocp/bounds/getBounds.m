function [bounds, scaling] = getBounds(Qs, GRF, model)

% This function assign the bounds to: model coordinates (Qs), Ground
% Reaction Forces (GRF), muscle activation, rate of change of muscle
% activation over time, MT force and rate of change of MT force over time. 
% Then, the bounds are scaled in such a way that they represents the +- 100%

%
% Inputs:
%   - Qs: 
%   - GRF:
%   - model:

import org.opensim.modeling.*

coordinateSet = model.getCoordinateSet();
muscleSet = model.getForceSet().getMuscles();
actuatorSet = model.getForceSet().getActuators();

num_q = coordinateSet.getSize();  % number of coordinates
num_muscles = muscleSet.getSize();  % number of muscles
num_act = actuatorSet.getSize();  % number of actuators

time = Qs.allfilt(:, 1);
y = zeros(num_q);

% Approximate 1st and 2nd derivative of Qs using analytical cubic spline
% derivation.
% loop through the columns of Qs, starting from index 2 (1 is time), and
% calculate the spline coefficients, given the experimental data
for i = 2:num_q
    cs = spline(time, Qs.allfilt(:, i));
    y(i) = eval_spline(cs, time, 2);
end

% ----------------------------- Qs bounds ----------------------------- %
% for the bounds, we follow the same index system of the coordinateSet
% object, which reflects the order of the coordinates that were assigned to
% the .osim model. We can just loop through the coordinateSet object.

% --- positions --- %
for i = 1:num_q
    bounds.Qs.upper(i) = max(y(i).pos);
    bounds.Qs.lower(i) = min(y(i).pos);
end
% scale
scaling.Qs = max(abs(bounds.Qs.lower),abs(bounds.Qs.upper));
bounds.Qs.lower = (bounds.Qs.lower)./scaling.Qs;
bounds.Qs.upper = (bounds.Qs.upper)./scaling.Qs;


% --- velocities --- %
for i = 1:num_q
    bounds.Qsdot.upper(i) = max(y(i).vel); 
    bounds.Qsdot.lower(i) = min(y(i).vel);
end
% scale
scaling.Qsdot = max(abs(bounds.Qsdot.lower),abs(bounds.Qsdot.upper));
bounds.Qsdot.lower = (bounds.Qsdot.lower)./scaling.Qsdot;
bounds.Qsdot.upper = (bounds.Qsdot.upper)./scaling.Qsdot;


% --- accelerations --- %
for i = 1:num_q
    bounds.Qsdotdot.upper(i) = max(y(i).acc);
    bounds.Qsdotdot.lower(i) = min(y(i).acc);
end
% scale
scaling.Qsdotdot = max(abs(bounds.Qsdotdot.lower),abs(bounds.Qsdotdot.upper));
bounds.Qsdotdot.lower = (bounds.Qsdotdot.lower)./scaling.Qsdotdot;
bounds.Qsdotdot.upper = (bounds.Qsdotdot.upper)./scaling.Qsdotdot;

% ----------------------------- GRFs bounds ----------------------------- %
lower_grf = min(GRF.val.all(:,2:end));
upper_grf = max(GRF.val.all(:,2:end));

% extend bounds to give some flexibility
force_range = abs(upper_grf - lower_grf);
bounds.GRF.lower = lower_grf - force_range;
bounds.GRF.upper = upper_grf + force_range;

% scale
scaling.GRF = max(abs(bounds.GRF.lower), abs(bounds.GRF.upper));
bounds.GRF.lower = (bounds.GRF.lower)./scaling.GRF;
bounds.GRF.upper = (bounds.GRF.upper)./scaling.GRF;


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