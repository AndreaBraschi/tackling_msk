function guess = xGuess(Qs, independent_coord_idx, time_mesh, time_col, num_q, scaling)
% --------------------------------------------------------------------------
% xGuess
%   This function sets the inital guess for X (q's and q_dot's
%   interleaved) from the experimental data.
%   
%   The experimental data is used to compute cubic spline coefficients, so
%   that a prediction on the values of q's can be made at the various
%   collocation points. Furthermore, velocity and acceleration of the
%   splines are computed at the same points to form a final initial guess.
%
%   This function was adapted to the specific case of this repo, but the
%   original, generic version can be found at:

%   https://github.com/KULeuvenNeuromechanics/PredSim/OCP/getGuess_DI_opti.m

% INPUTs:
%   - Qs (struct): q's coming IK.
%   - independent_coord_idx (double): array containing indices of Qs that correspond to the model independent coordinates.
%   - num_q (int): number of coordinates
%   - scaling (struct): object containing the scaling factors for all the
%   state vector dimensions.
 

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
    validate_array(y.pos, "Q_spline")

    Qdots_spline(:, i) = y.vel;
    validate_array(y.vel, "Qdots_spline")
    
    
    Qdotdots_spline(:, i) = y.acc;
    validate_array(y.acc, "Qdotdots_spline")


    
    % evaluate spline (and compute 1st and 2nd derivative) at each 
    % collocation points
    y = eval_spline_col(cs, time_col, col_k, col_dt, 2);

    Qs_spline_col(:, i) = y.pos;
    validate_array(y.pos, "Qs_spline_col")

    
    Qdots_spline_col(:, i) = y.vel;
    validate_array(y.vel, "Qdots_spline_col")
    
    Qdotdots_spline_col(:, i) = y.acc;
    validate_array(y.acc, "Qdotdots_spline_col")
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
guess.Qdotdots_col_ind = Qdotdots_spline_col(:, independent_coord_idx);


% Do the same for the independent coordinates only.
% We first need to place Qs and Qsdot as Simbody/OpenSim expect the state
% vector to be: Q = [q_dot(:, 1), q_dot(:, 1), q(:, 2), q_dot(:, 2), ...]
dims_ind = size(independent_coord_idx, 2) * 2;
Q_ind = cat(3, Qs_spline(:, independent_coord_idx), Qdots_spline(:, independent_coord_idx));  % [T_mesh x num_q x 2]
Q_ind = reshape(permute(Q_ind, [1, 3, 2]), T_mesh, dims_ind);

Q_col_ind = cat(3, Qs_spline_col(:, independent_coord_idx), Qdots_spline_col(:, independent_coord_idx));  % [T_col x num_q x 2]
Q_col_ind = reshape(permute(Q_col_ind, [1, 3, 2]), T_col, dims_ind);

% add to a 'guess' struct: we can add the Q acceleration as they are, as
% acceleration isn't part of the state vector.

% end points of the mesh
guess.Qs_all_ind = Q_ind;

% collocation points
guess.Qs_col_ind = Q_col_ind;


% look for invalid numerical values
validate_guess(guess);


end