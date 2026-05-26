function guess = getGuess(Qs, time_mesh, time_col, num_q, num_muscles, num_actuators, scaling)
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
%   - num_muscles (int): number of muscles
% 
%   - num_actuators (int): number of actuators
%
%   - scaling (struct): object containing the scaling factors for all the
%   state vector dimensions.
%
% OUTPUT:
%   - guess (struct): guess object that contains the initial guess for all
%   the dimension of the state vector.
 

T_mesh = size(time_mesh, 2);
T_col = size(time_col, 1);
dims = num_q * 2;  % total dimensionality of the q-part of the state vector.

% --------------------------------------------------------------------------
% ------ Interpolate the generalised coordinates using cubic splines ------
% Pre allocate arrays for
% Position
Qs_spline = zeros(size(Qs.allinterpfilt));
Qs_spline_col = zeros(size(Qs.allinterpfilt_col));

% Velocity
Qdots_spline = zeros(size(Qs.allinterpfilt));
Qdots_spline_col = zeros(size(Qs.allinterpfilt_col));

% Acceleration
Qdotdots_spline = zeros(size(Qs.allinterpfilt));
Qdotdots_spline_col = zeros(size(Qs.allinterpfilt_col));

mesh_k = discretize(time_mesh', Qs.time);
mesh_dt = time_mesh' - Qs.time(mesh_k);

col_k = discretize(time_col, Qs.time);
col_dt = time_col - Qs.time(col_k);

for i = 1:num_q
    % calculate T - 1 (T being period) spline coefficients for Qs
    cs = spline(Qs.time, Qs.allfilt(:, i + 1));  % i + 1, because i = 1 is time column
    
    % evaluate spline (and compute 1st and 2nd derivative) at the end points 
    % of each mesh segment
    y = eval_spline_col(cs, time_mesh, mesh_k, mesh_dt, 2);

    Qs_spline(:, i) = y.pos;
    Qdots_spline(:, i) = y.vel;
    Qdotdots_spline(:, i) = y.acc;


    
    % evaluate spline (and compute 1st and 2nd derivative) at each 
    % collocation points
    y = eval_spline_col(cs, time_col, col_k, col_dt, 2);

    Qs_spline_col(:, i) = y.pos;
    Qdots_spline_col(:, i) = y.vel;
    Qdotdots_spline_col(:, i) = y.acc;

end


% ----- scale ----- %
% end-points
Qs_spline = (Qs_spline)./scaling.Qs;
Qdots_spline = (Qdots_spline)./scaling.Qsdot;
Qdotdots_spline = (Qdotdots_spline)./scaling.Qsdotdot;

% collocation points
Qs_spline_col = (Qs_spline_col)./scaling.Qs;
Qdots_spline_col = (Qdots_spline_col)./scaling.Qsdot;
Qdotdots_spline_col = (Qdotdots_spline_col)./scaling.Qsdotdot;

% add the splined Q, Qdot and Qdotdot to a 'guess' struct

% We first need to place Qs and Qsdot as Simbody/OpenSim expect the state
% vector to be: Q = [q_dot(:, 1), q_dot(:, 1), q(:, 2), q_dot(:, 2), ...]
Q = cat(3, Qs_spline, Qdots_spline);  % [T_mesh x num_q x 2]
Q = reshape(permute(Q, [1, 3, 2]), T_mesh, dims);

Q_col = cat(3, Qs_spline_col, Qdots_spline_col);  % [T_col x num_q x 2]
Q_col = reshape(permute(Q_col, [1, 3, 2]), T_col, dims);

% add to a 'guess' struct: we can add the Q acceleration as they are, as
% acceleration isn't part of the state vector.

% end points of the mesh
guess.Qs_all = Q;
guess.Qdotdots_all = Qdotdots_spline;

% collocation points
guess.Qs_col = Q_col;
guess.Qdotdots_col = Qdotdots_spline_col;


% ----- Muscle variables ----- %
% mesh end-points
a = 0.1 * ones(T_mesh, num_muscles);
vA = 0.01 * ones(T_mesh, num_muscles);
FTtilde = 0.1 * ones(T_mesh, num_muscles);
dFTtilde = 0.01 * ones(T_mesh, num_muscles);

% scale
guess.a = (a)./scaling.a;
guess.FTtilde = FTtilde./scaling.FTtilde;
guess.vA  = (vA)./scaling.vA;
guess.dFTtilde = (dFTtilde)./scaling.dFTtilde;

% collocation points
a_col = 0.1 * ones(T_col, num_muscles);
FTtilde_col = 0.1 * ones(T_col, num_muscles);
dFTtilde_col = 0.01 * ones(T_col, num_muscles);

% scale
guess.a_col = (a_col)./scaling.a;
guess.FTtilde_col = FTtilde_col./scaling.FTtilde;
guess.dFTtilde_col = (dFTtilde_col)./scaling.dFTtilde;


% Torque actuators
guess.a_a = 0.1 * ones(T_mesh, num_actuators);
guess.e_a = 0.1 * ones(T_col, num_actuators);
guess.a_a_col= 0.01 * ones(T_mesh, num_muscles);



end