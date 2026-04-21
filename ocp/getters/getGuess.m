function guess = getGuess(Qs, nq, nMuscles, nActuators)
% --------------------------------------------------------------------------
% getGuess
%   This function creates an inital guess for the design variables directly
%   from the experimental data. The design variable correspond to the
%   dimensions of the state vector, which is made of: q, q_dot, muscle act,
%   muscle dact/dt, 
%   This function was adapted to the specific case of this repo, but the
%   original, generic version can be found at:

%   https://github.com/KULeuvenNeuromechanics/PredSim/OCP/getGuess_DI_opti.m

% INPUTs:
%   - Qs (struct): coordinates coming from IK.
% 
%   - nq (int): number of coordinates
%
%   - nMuscles (int): number of muscles
% 
%   - nActuatos (int): number of actuators
%
% OUTPUT:
%   - guess (struct): guess object that contains the initial guess for all
%   the dimension of the state vector.
 

T = size(Qs.allinterpfilt, 1);  % Period: discrete points along the trajectory
dims = nq * 2;  % total dimensionality of the q-part of the state vector.

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

% We first need to place Qs and Qsdot as Simbody/OpenSim expect the state
% vector to be: Q = [q_dot(:, 1), q_dot(:, 1), q(:, 2), q_dot(:, 2), ...]
Q = reshape([Qs_spline, Qdots_spline], T, 2, dims);
Q = reshape(permute(Q, [1, 3, 2]), T, dims * 2);

Q_col = reshape([Qs_spline_col, Qdots_spline_col], T, 2, dims);
Q_col = reshape(permute(Q_col, [1, 3, 2]), T, dims * 2);

% add to a 'guess' struct: we can add the Q acceleration as they are, as
% acceleration isn't part of the state vector.

% end points of the mesh
guess.Qs_all = Q;
guess.Qdotdots_all = Qdotdots_spline.data;

% collocation points
guess.Qs_col = Q_col;
guess.Qdotdots_col = Qdotdots_spline_col.data;


% Muscle variables
guess.a = 0.1 * ones(T, NMuscle);
guess.vA = 0.01 * ones(T, NMuscle);
guess.FTtilde = 0.1 * ones(T, NMuscle);
guess.dFTtilde = 0.01 * ones(T, NMuscle);


% Torque actuators
guess.a_a = 0.1*ones(T, nActuators);
guess.e_a = 0.1*ones(T, nActuators);



end