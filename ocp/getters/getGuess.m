function guess = getGuess(Qs, nq, nMuscles, nActuators)
% --------------------------------------------------------------------------
% getGuess
%   This function creates an inital guess for the design variables directly
%   from the experimental data. The design variable correspond to the
%   dimensions of the state vector, which is made of: q, q_dot, muscle act,
%   muscle dact/dt, 
%   This function was adapted to the specific case of this repo, but the
%   original form can be found at:

%   https://github.com/KULeuvenNeuromechanics/PredSim/OCP/getGuess_DI_opti.m

% INPUTs:
%   - Qs: coordinates coming from IK.
% 
%   - nq (integer): number of coordinates

%   - scaling -
% 
%   - d: degree of the interpolating polynomial of the collocation scheme
%
% OUTPUT:
%   - guess -
%   * initial guess values for all optimisation variables
% 

T = size(Qs.allinterpfilt, 1);  % Period: discrete points along the trajectory


% --------------------------------------------------------------------------
% ------ Interpolate the generalised coordinates using cubic splines ------
% Pre allocate arrays for
% Position
Qs_spline.data = zeros(size(Qs.allinterpfilt));
Qs_spline.data(:,1) = Qs.allinterpfilt(1:end, 1);

% Velocity
Qdots_spline.data = zeros(size(Qs.allinterpfilt));
Qdots_spline.data(:,1) = Qs.allinterpfilt(1:end,1);

% Acceleration
Qdotdots_spline.data = zeros(size(Qs.allinterpfilt));
Qdotdots_spline.data(:,1) = Qs.allinterpfilt(1:end,1);


for i = 2:size(Qs.allfilt,2)
    % calculate T - 1 (T being period) spline coefficients for Qs
    Qs.datafiltspline(i) = spline(Qs.allfilt(:,1), Qs.allfilt(:,i));
    
    % evaluate spline (and compute 1st and 2nd derivative) at the end points 
    % of each mesh segment
    [Qs_spline.data(:,i), Qdots_spline.data(:,i), Qdotdots_spline.data(:,i)] = ...
        eval_spline(Qs.datafiltspline(i), Qs.allinterpfilt(1:end,1),1);
    
    % evaluate spline (and compute 1st and 2nd derivative) at each 
    % collocation points
    [Qs_spline_col.data(:,i), Qdots_spline_col.data(:,i), Qdotdots_spline_col.data(:,i)] = ...
        SplineEval_ppuval(Qs.datafiltspline(i), Qs.allinterpfilt_col(1:end,1), 1);
end


% add the splined Q, Qdot and Qdotdot to a 'guess' struct
for i=1:nq
    % end points of the mesh
    guess.Qs_all(:, i) = Qs_spline.data(:, i);
    guess.Qdots_all(:, i) = Qdots_spline.data(:, i);
    guess.Qdotdots_all(:, i) = Qdotdots_spline.data(:, i);

    % collocation points
    guess.Qs_col(:, i) = Qs_spline_col.data(:, i);
    guess.Qdots_col(:, i) = Qdots_spline_col.data(:, i);
    guess.Qdotdots_col(:, i) = Qdotdots_spline_col.data(:, i);
end


% Muscle variables
guess.a = 0.1 * ones(T, NMuscle);
guess.vA = 0.01 * ones(T, NMuscle);
guess.FTtilde = 0.1 * ones(T, NMuscle);
guess.dFTtilde = 0.01 * ones(T, NMuscle);


% Torque actuators
guess.a_a = 0.1*ones(T, nActuators);
guess.e_a = 0.1*ones(T, nActuators);



end