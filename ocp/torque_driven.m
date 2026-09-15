function [opti, w_opt, stats, g_opt, lambda_x, lambda_g] = torque_driven( ...
    model_path, ik_path, dll_path, config_filepath, kinematic_coupling_path, ...
    output_dir, grf_path ... 
    )

% -------------------------------------------------------------------------
% torque_driven
%   This function aims at tracking experimental data while solving an
%   Optimal Control Problem using idealised torque actuator only (no muscles).

% INPUTs:
%   - model_path (str): path to the OpenSim model (.osim file)

%   - ik_path (str): path to the IK solution (.mot file) 

%   - grf_path (str): path to the Ground Reaction Force file (optional).

%   - dll_filepath (char): path to the dll specific to the current trial 
%     being tracked. char MUST be the data type expected by CasADi.
% -------------------------------------------------------------------------
    arguments
        model_path
        ik_path
        dll_path
        config_filepath
        kinematic_coupling_path
        output_dir
        grf_path = []   % default to empty if not passed
    end        

    import casadi.*
    import org.opensim.modeling.*
    
    
    % read project config file
    json_str = fileread(config_filepath);
    config_struct = jsondecode(json_str);

   
    % parallelisation settings
    parallelMode = 'thread';
    num_threads = 8; % Number of threads used in parallel.

    
    % --- add paths of subdirectories where we'll be picking functions from --- %     
    % various utility functions
    pathUtils = [pwd,'/utils'];
    addpath(genpath(pathUtils));

    % collocation scheme
    pathCollocationScheme = [pwd,'/collocationScheme'];
    addpath(genpath(pathCollocationScheme));

    pathCasADiFunctions = [pwd,'/casadi_functions'];
    addpath(genpath(pathCasADiFunctions));

    pathBounds = [pwd,'/bounds'];
    addpath(genpath(pathBounds));

    pathGuess = [pwd,'/guess'];
    addpath(genpath(pathGuess));

    pathGetters = [pwd,'/getters'];
    addpath(genpath(pathGetters));

    pathConstraints = [pwd,'/kinematic_coupling'];
    addpath(genpath(pathConstraints));

    pathOpt = [pwd,'/opt'];
    addpath(genpath(pathOpt));
    
    pathEval = [pwd,'/evaluate'];
    addpath(genpath(pathEval));

    pathMuscleModel = [pwd,'/muscle_model'];
    addpath(genpath(pathMuscleModel));   
    
   

    % Collocation scheme
    N = config_struct.("collocation").("number_of_segments");   % number of mesh intervals
    d = config_struct.("collocation").("num_points"); % number of collocation points per mesh interval
    method = config_struct.("collocation").("method"); % collocation method
    
    [tau_root, C, D, B] = collocationScheme(d, method); % collocation scheme.
    
    
    % Load external functions
    F = external('F', dll_path); 
    % load information relevant to external function from config file
    dll_force_indices = config_struct.("dll_force_indices");
    dll_grf_indices = [dll_force_indices.rGRF', dll_force_indices.lGRF'];

    % -------------------------- OpenSim Model -------------------------- % 
    model = Model(model_path);

     % ---------- sets ----------- %
    coordinateSet = model.getCoordinateSet();

    % ---------- info ----------- %
    % get coordinate names - all: independent and dependent
    q_names = getItemNames(coordinateSet);
    num_q_all = size(q_names, 2);
    q_all_indices = 1:num_q_all;

    % read names of dependent coordinates
    q_dep_names = config_struct.("dependent_coord_names");
    
    % find dependent coordinates indices
    q_dep_indices = cellfun(@(name) find(strcmp(q_names, name)), q_dep_names);

    q_dep_names_check = q_names(q_dep_indices);

    % independent coordinates indices
    q_ind_indices = setdiff(q_all_indices, q_dep_indices);
    q_ind_names = q_names(q_ind_indices);

    % differentiate between number of independent and dependent coords 
    num_q_dep = size(q_dep_names, 1);
    num_q_ind = num_q_all - num_q_dep;   

    % check if the user has defined some special name not to be tracked
    no_tracking_names = config_struct.tracking.no_tracking;
    special_tracking_names = config_struct.tracking.special_tracking;

    % check if the user has set some special bounds to be assigned
    if ~isempty(no_tracking_names)
        no_tracking_indices = cellfun(@(name) find(strcmp(q_names, name)), no_tracking_names);
        num_no_tracking_indices = size(no_tracking_names, 1);

        q_tracking_indices = setdiff(q_ind_indices, no_tracking_indices');
        
        % now, let's find the indices of the special DoFs within
        % the set of independent DoFs
        no_tracking_indices_ind = cellfun(@(name) find(strcmp(q_ind_names, name)), no_tracking_names');
        local_ind_indices = 1:num_q_ind;
        local_tracking_indices = setdiff(local_ind_indices, no_tracking_indices_ind);

    else
        local_ind_indices = 1:num_q_ind;
        local_tracking_indices = local_ind_indices;
        q_tracking_indices = q_ind_indices;

    end

    num_q_tracking = size(q_tracking_indices, 2);
    tracking_indices = 1:num_q_tracking;
    
    % check if the user has decleared a specific q to be tracked
    % differently (i.e. different weight: higher or lower)
    if ~isempty(special_tracking_names)
        special_indices = cellfun(@(name) find(strcmp(q_names, name)), special_tracking_names);
        num_special_indices = size(special_indices, 1);
        
        % now, let's find the indices of the special DoFs within
        % the set of independent DoFs
        special_indices_ind = cellfun(@(name) find(strcmp(q_ind_names, name)), special_tracking_names');
        indices_to_track = setdiff(tracking_indices, special_indices_ind);
    else
        indices_to_track = tracking_indices;

    end

    
    
    % ------------------------ Experimental Data ------------------------ %
    % import scaled factors
    scale_factors = config_struct.actuator_scale_factors;
    joint_moments_sf = scale_factors.("joint_moments");
    
    res_sf = scale_factors.residuals;
    res_keys      = fieldnames(res_sf);
    res_indices = [];
    for i = 1:length(res_keys)
        key   = res_keys{i};
        items = res_sf.(key);        

        res_sf_value.(key) = items.magnitude;
        res_sf_names = items.names;
        res_sf_indices.(key) = cellfun(@(name) find(strcmp(q_names, name)), res_sf_names);
        res_indices = [res_indices; res_sf_indices.(key)];
    end

    % pelvis_y_sf = scale_factors.("pelvis_y");

    
    % ------ IK ------ %
    Qs = readMotQs(ik_path, 15);
    % Qs.allfilt(:, pelvis_y_idx + 1) = Qs.allfilt(:, pelvis_y_idx + 1) + ...
    %     Qs.allfilt(:, pelvis_y_idx + 1) * pelvis_y_sf;


    
    % zero the values of Qs corresponding to the spine translational DoFs.
    % We can do it here as we're not going to track the experimental data,
    % and we're going to make our life easier when it comes to guess
    % generation (0) and bounds ([-1 1])
    if ~isempty(no_tracking_names)
        Qs.allfilt(:, no_tracking_indices + 1) = 0;

    end
    
 
    % we define one actuator per independent DoF
    % here we check which coordinates the user has indicated as where to
    % assign the residual forces. We need to make sure that we correctly
    % separate the residual forces from the actuator forces/torques so that
    % we can assign them to the correct set of constraints.
    actuator_indices = setdiff(q_ind_indices, res_indices);
    if ~isempty(no_tracking_names)
        actuator_indices = setdiff(actuator_indices, no_tracking_indices);

    end
    % actuator_indices = actuator_indices(:, 1:5);
    actuator_names = q_names(actuator_indices);
    num_actuators = size(actuator_indices, 2);

    % prepare actuator scale factor vector to be used later in the equality
    % constraint set
    actuator_sf = ones(num_actuators, 1) * joint_moments_sf;

    % let's check whether the user has provided some actuators to be given
    % a specific scale factor
    act_special_sf = scale_factors.special;
    act_keys      = fieldnames(act_special_sf);
    for i = 1:length(act_keys)
        key   = act_keys{i};
        items = act_special_sf.(key);        

        act_special_value = items.magnitude;
        act_special_names = items.names;
        
        act_special_indices = cellfun(@(name) find(strcmp(actuator_names, name)), act_special_names);
        actuator_sf(act_special_indices, :) = act_special_value;
    end


    % see if the user provided any GRF
    if ~isempty(grf_path)
        GRFs = readMotGrf(grf_path, 15);
        experimental_force_indices = config_struct.("experimental_force_indices");
        grf_indices = [experimental_force_indices.rGRF', experimental_force_indices.lGRF'];
    end



    % read initial and final time from IK 
    time_opt = [Qs.time(1, 1) Qs.time(end, 1)];


    % -------------------------- Interpolation -------------------------- %
    % compute the Period of each mesh
    mesh_T = (time_opt(2) - time_opt(1)) / N;

    % create a time vector that represents N + 1 points which are the last
    % point of each mesh. Considering the continuity constraint, the last
    % point of one mesh represents the first of the next one.
    time_intervals = time_opt(1):mesh_T:time_opt(2);
    
    % now, for each mesh, let's create a vector that contains the 3
    % collocation points along the input space.
    time_grid = zeros(N * d, 1);
    time_grid(1:d:end) = time_intervals(1:end-1) + tau_root(2) * mesh_T; % 1st collocation point of every mesh
    time_grid(2:d:end) = time_intervals(1:end-1) + tau_root(3) * mesh_T; % 2nd collocation point of every mesh
    time_grid(3:d:end) = time_intervals(1:end-1) + tau_root(4) * mesh_T; % 3rd collocation point of every mesh


    % --- IK --- %
    % find Qs values at first/last point of each mesh
    Qs.allinterpfilt = interp1(Qs.time, Qs.allfilt(:, 2:end), time_intervals);
    
    % find Qs values at each collocation point along the trajectory
    Qs.allinterpfilt_col = interp1(Qs.time, Qs.allfilt(:, 2:end), time_grid');

    
    % ----------------------------- Bounds  ----------------------------- %
    % ---------- X ----------- %
    special_qs = config_struct.bounds.special_qs;
    [x_bounds, x_scaling] = set_xBounds(Qs, q_ind_indices, num_q_all, special_qs);

    % ---------- actuators ----------- %
    [act_bounds, act_scaling] = set_torqueBounds(num_actuators);
    
    % ---------- residuals ----------- %
    residual_bounds = set_residualBounds(config_struct);

    % repeat bound values over time points
    lb_a_a = ones(N + 1, 1) * act_bounds.a_a.lower;
    ub_a_a = ones(N + 1, 1) * act_bounds.a_a.upper;

    lb_a_a_col = ones(N * d, 1) * act_bounds.a_a.lower;
    ub_a_a_col = ones(N * d, 1) * act_bounds.a_a.upper;
        
    % -------------------------- Initial Guess  ------------------------- %
    torque_guess = torqueGuess(num_actuators, time_intervals, time_grid);
    x_guess = xGuess(Qs, q_ind_indices, time_intervals, time_grid, ...
            num_q_all, x_scaling);
    

    save(fullfile(output_dir, "guess.mat"), "x_guess", "torque_guess");

    
    % check whether initial guess are outside of the bounds
    qs_guess = x_guess.Qs_all_ind;
    qs_lower = x_bounds.X.lower;
    qs_upper = x_bounds.X.upper;
    interleaved = reshape(repmat(q_ind_names, 2, 1), 1, []);
    for n = 1:size(qs_guess, 2)
        q_guess = qs_guess(:, n);
        q_lower = qs_lower(:, n);
        q_upper = qs_upper(:, n);
        var_name = interleaved{n}; 

        if any(q_guess > q_upper)
            excess = max(q_guess - q_upper);
            fprintf('Variable %s exceeds upper bound at n = %d by %.4f.\n', var_name, n, excess);

        end
        if any(q_guess < q_lower)
            excess = max(q_lower - q_guess);
            fprintf('Variable %s exceeds lower bound at n = %d by %.4f.\n', var_name, n, excess);
        end

    end


    qs_guess = x_guess.Qdotdots_col_ind;
    qs_lower = x_bounds.Qsdotdot.lower_ind;
    qs_upper = x_bounds.Qsdotdot.upper_ind;
    for n = 1:size(qs_guess, 2)
        q_guess = qs_guess(:, n);
        q_lower = qs_lower(:, n);
        q_upper = qs_upper(:, n);
        var_name = q_names{n}; 

        if any(q_guess > q_upper)
            excess = max(q_guess - q_upper);
            fprintf('Variable %s exceeds upper bound at n = %d by %.4f.\n', var_name, n, excess);

        end
        if any(q_guess < q_lower)
            excess = max(q_lower - q_guess);
            fprintf('Variable %s exceeds lower bound at n = %d by %.4f.\n', var_name, n, excess);
        end

    end
    
    
    % ------------------- Experimental Data Scaling  ------------------- %
    % we have the scaled experimental data stored in the 'guess' struct.
    % q
    Qs_scaled_all = x_guess.Qs_all(:, 1:2:end);
    Qs_scaled_col_all = x_guess.Qs_col(:, 1:2:end);
        
    % retrieve just independent coordinates
    Qs_scaled = Qs_scaled_all(:, q_tracking_indices);
    Qs_scaled_col = Qs_scaled_col_all(:, q_tracking_indices);
    
    % q dot
    Qdots_scaled_all = x_guess.Qs_all(:, 2:2:end);
    Qdots_scaled_col_all = x_guess.Qs_col(:, 2:2:end);
    % retrieve just independent coordinates
    Qdots_scaled = Qdots_scaled_all(:, q_tracking_indices);
    Qdots_scaled_col = Qdots_scaled_col_all(:, q_tracking_indices);

    save(fullfile(output_dir, "x_scaled_col.mat"), "Qs_scaled_col", "Qdots_scaled_col", "time_grid", "mesh_T");
    

    save(fullfile(output_dir, "scaling.mat"), "x_scaling", "act_scaling");

    config_out = struct();
    config_out.q_names              = q_names;
    config_out.q_ind_indices        = q_ind_indices;
    config_out.num_q_ind            = num_q_ind;
    config_out.num_actuators        = num_actuators;
    config_out.d                    = d;
    config_out.N                    = N;
    config_out.num_q_all            = num_q_all;
    config_out.local_tracking_indices           = local_tracking_indices;
    config_out.indices_to_track     = indices_to_track ;
    config_out.actuator_indices     = actuator_indices ;
    config_out.actuator_sf          = actuator_sf ;
    config_out.num_res_vars          = size(res_indices, 1) ;
    
    if ~isempty(special_tracking_names)  
        config_out.special_indices_ind  = special_indices_ind;
    end

    % term for special DoFs
    if ~isempty(no_tracking_names)
        config_out.no_tracking_indices_ind  = no_tracking_indices_ind;
    end
    
    
    json_str = jsonencode(config_out, 'PrettyPrint', true);
    fid = fopen(fullfile(output_dir, 'opt_setup.json'), 'w');
    fprintf(fid, '%s', json_str);
    fclose(fid);


    % ------------------------------------------------------------------- %
    %                          NLP formulation                            %
    % ------------------------------------------------------------------- %

    % Start with an empty NLP. Initialize opti instance.
    opti = casadi.Opti();

    % Note: we're trying to optimise the state variable and the controls at
    % the:
    % 1) collocation points
    % 2) mesh end-points: this is necessary to enforce continuity across
    % segments
    % We follow the same scheme for every subset of the design variable.

    % ----- Qs and Qd_dot ----- %
    % 1) collocation points
    % Create a symbolic variable within the problem and assign it a number
    % of dimensions and points along the trajectory.
    dims = 2 * num_q_ind;
    points = N * d;  % we have d-collocation points x N-segments
    X_col = opti.variable(dims, points);
    
    % bounds
    opti.subject_to(x_bounds.X.lower' < X_col < x_bounds.X.upper');
    % initial condition
    opti.set_initial(X_col, x_guess.Qs_col_ind');
    
    % 2) mesh end-points
    points = N + 1;  
    X = opti.variable(dims, points);
    opti.subject_to(x_bounds.X.lower' < X < x_bounds.X.upper');
    opti.set_initial(X, x_guess.Qs_all_ind');


    % ----- Torque Actuators ----- %
    dims = num_actuators;
    points = d * N;

    % 1) collocation points
    a_a_col = opti.variable(dims, points);
    opti.subject_to(lb_a_a_col'< a_a_col < ub_a_a_col');
    opti.set_initial(a_a_col, torque_guess.a_a_col');

    % 2) mesh end-points
    points = N + 1;
    a_a = opti.variable(dims, points);
    opti.subject_to(lb_a_a'< a_a < ub_a_a');
    opti.set_initial(a_a, torque_guess.a_a');  


    fprintf('number of states   : %d\n', num_q_ind * 2 + num_actuators * 2);

    % ----------------------- Controls  ----------------------- %    
    % ----- Actuator Excitation ----- %
    dims = num_actuators;
    e_a = opti.variable(dims, points);
    opti.subject_to(act_bounds.e_a.lower' < e_a < act_bounds.e_a.upper');
    opti.set_initial(e_a, torque_guess.e_a');
    % 

    
    % Time derivative of Qdots (states) at collocation points
    A_col = opti.variable(num_q_ind, d * N);
    opti.subject_to(x_bounds.Qsdotdot.lower_ind' < A_col < x_bounds.Qsdotdot.upper_ind');
    opti.set_initial(A_col, x_guess.Qdotdots_col_ind'); 

    fprintf('number of controls : %d\n', num_actuators + num_q_ind);


    func_map_in = {X(:, 1:end-1), X_col, A_col, ...
        a_a(:, 1:end-1), a_a_col, e_a(:, 1:end-1), ...
        MX(Qs_scaled(1:end-1, :)'), MX(Qs_scaled_col'), ...
        MX(Qdots_scaled(1:end-1, :)'), MX(Qdots_scaled_col')};
    

    if ~isempty(grf_path)
        % --- GRF --- %
        % find indices where the 2 items of time_opt are 
        dt_GRF = GRFs.time(2) - GRFs.time(1);
        grf_init = find((GRFs.time<(time_opt(1) + dt_GRF/2)) & (GRFs.time>=(time_opt(1) - dt_GRF/2)));
        grf_end = find((GRFs.time<(time_opt(2) + dt_GRF/2)) & (GRFs.time>=(time_opt(2) - dt_GRF/2)));
    
        % crop force data
        GRFs.time = GRFs.time(grf_init:grf_end);
        GRFs.data = GRFs.data(grf_init:grf_end, :);
    
        % find values at:
        % collocation points
        GRF_col = interp1(GRFs.time, GRFs.data(:, grf_indices), time_grid');

        save(fullfile(output_dir, "GRF_col.mat"), "GRF_col");

        scaling_GRF = max(abs(min(GRF_col)), abs(max(GRF_col)));
        GRF_scaled = GRF_col./scaling_GRF;

        func_map_in{end+1} = MX(GRF_scaled');

    end
    
    
    
    % ----------------------- Residuals  ----------------------- %
    % we use the config file to register the residuals that the user
    % specified:
    residuals = config_struct.bounds.residuals;
    res_keys      = fieldnames(residuals);

    for i = 1:length(res_keys)
        key   = res_keys{i};
        values = residuals.(key);        
        var_size = values(2);       % how many DoFs it is applied to
        
        % register design variable
        residual_vars.(key) = opti.variable(var_size, d * N);
        opti.subject_to(residual_bounds.(key).lower'< residual_vars.(key) < residual_bounds.(key).upper');
        opti.set_initial(residual_vars.(key), zeros(var_size, N * d));

        func_map_in{end+1} = residual_vars.(key);
            
    end
    

    % print statement to check whether the number of design variables
    % matches up with the number of bounds.
    fprintf('Decision variables : %d\n', opti.nx);
    fprintf('Bound Constraints        : %d\n', opti.ng);


    % ------------------------------------------------------------------- %
    % The following section uses the "MX" object from the CasADi library
    % to register symbolic variables that will later be used in CasADi 
    % "Function" objects, which dictate how these variables interact with
    % each other.

    % The suffix "k" indicates the variable at the mesh end-points.
    % "j", on the other hand, at the collocation points.
    % ------------------------------------------------------------------- %

    % ---------- state-vector ---------- %   
    % q and q_dot
    Xk = MX.sym('Xk', 2 * num_q_ind);   % shape: [num_q * 2]
    Xj = MX.sym('Xj', 2 * num_q_ind, d); % shape: [num_q * 2, d]
    Xkj = [Xk Xj];  % shape: [num_q * 2, d + 1]
    
    % torque values
    a_ak = MX.sym('a_ak', num_actuators);
    a_aj = MX.sym('a_aj', num_actuators, d);
    a_akj = [a_ak a_aj];

    
    % ---------- controls ---------- %

    % torque actuator excitation
    e_ak = MX.sym('e_ak', num_actuators);

    
    % q_dotdot: remember, we're using implicit formulation. Therefore,
    % accelerations are treated as "controls".
    Aj = MX.sym('Aj', num_q_ind, d);

    % ------- Experimental Data to Track ------- %
    % ------------------------ TODO -----------------------------------%
    % here we need to update the number of q's that is used to initialise
    % the experimental data to track
    Qs_track_k = MX.sym('Qs_track_k', num_q_tracking);
    Qs_track_j = MX.sym('Qs_track_j', num_q_tracking, d);
    Qs_track_kj = [Qs_track_k Qs_track_j];
    
    Qdots_track_k = MX.sym('Qdots_track_k', num_q_tracking);
    Qdots_track_j = MX.sym('Qdots_track_j', num_q_tracking, d);
    Qdots_track_kj = [Qdots_track_k Qdots_track_j];

    func_in = {Xk, Xj, Aj, a_ak, a_aj, e_ak, Qs_track_k, Qs_track_j, ...
        Qdots_track_k, Qdots_track_j};
    
   
    if ~isempty(grf_path)
        num_grfs = size(GRF_scaled, 2);
        GRF_track_j = MX.sym('GRF_track_j', num_grfs, d); 

        func_in{end+1} = GRF_track_j;

    end
    
    
    % % ---------- Residuals ---------- %
    for i = 1:length(res_keys)
        key   = res_keys{i};
        values = residuals.(key);        
        var_size = values(2);   
        
        res_sym.(key) = MX.sym([key, '_j'], var_size, d);

        func_in{end+1} = res_sym.(key);
            
    end
    % initialise set of constraint vector
    eq_constr = {}; % equality constraint vector

    % ------------------------------------------------------------------- %
    % The following section is where the external function is called.
    % Furthermore, here is where the CasADi function create the
    % symbolic interdependency between the variables that were previously
    % registered.

    % We do it in such a way that can be parallelised across all the
    % trajectory segments. This means that the trajectory segments are
    % treated independently from one another and that the following
    % operation can be parallelised across multiple threads.
    % ------------------------------------------------------------------- %
    % we need to make use of the unscaled version of the data
    Xkj_nsc = MX.zeros(size(Xkj));
    Xkj_nsc(1:2:end, :) = Xkj(1:2:end,:).*x_scaling.Qs(:, q_ind_indices)';
    Xkj_nsc(2:2:end, :) = Xkj(2:2:end,:).*x_scaling.Qsdot(:, q_ind_indices)';
    
    Aj_nsc = Aj.*(x_scaling.Qsdotdot(:, q_ind_indices)');  
    
  
    % ---------- CasADi functions ---------- %
    % Torque actuation dynamics
    activation_dynamics_function = torque_activation_dynamics_casadi(pathMuscleModel, num_actuators);

    % --- sum of squares --- %
    % q
    J_q = sum_of_squares('q', num_q_tracking - num_special_indices);
    
    % q_dot
    J_q_dot = sum_of_squares('q_dot', num_q_tracking);

    % q_dot_dot
    J_acc = sum_of_squares('acc', num_q_ind);

    % special q-tracking term that the user has defined
    % we initialise this anyway so that even if the user hasn't provided
    % any, we just add a 0 to the objective.
    q_no_tracking_term = MX.zeros(1, 1);
    if ~isempty(no_tracking_names)
        J_no_tracking = sum_of_squares('no_tracking', num_no_tracking_indices);

    end

    q_special_term = MX.zeros(1, 1);
    if ~isempty(special_tracking_names)
        J_special = sum_of_squares('special', num_special_indices);

    end

    GRF_term = MX.zeros(1, 1);
    if ~isempty(grf_path)
        J_GRF = sum_of_squares('GRF', num_grfs);
    end

    
    
    for i = 1:length(res_keys)
        key   = res_keys{i};
        values = residuals.(key);        
        var_size = values(2);       % how many DoFs it is applied to
        
        % register residual cost functions
        J_res.(key) = sum_of_squares(key, var_size);
        res_term.(key) = MX.zeros(1, 1);
           
    end
    
    fprintf("\nAll CasADi functions have been registered\n");

    
    % Read the weights of the cost function from config file
    W = config_struct.("W");
    
    % ------------------------------------------------------------------- %
    % Preparae accumulators for what we want to numerically evaluate after
    % optimisation:

    % cost function
    J = MX.zeros(1, 1);
    q_term = MX.zeros(1, 1);
    q_dot_term = MX.zeros(1, 1);
    acc_term = MX.zeros(1, 1);
    residuals_term = MX.zeros(1, 1);

    
    fprintf('\nCollocation loop starting\n')
    % for one segment:
    % loop through the collocation points
    for i = 1:d 
        
        % current q-part of state vector
        x_i = Xkj_nsc(:, i + 1);

        % current q accelerations
        acc_i = Aj_nsc(:, i);

        [x_all_i, acc_all_i] = apply_constraints(x_i, acc_i, q_names, q_ind_indices, kinematic_coupling_path);

        % evaluate external function
        Ti = F(vertcat(x_all_i, acc_all_i));
       

        % ---------------- Cost Function Terms ---------------- %
        % Q
        % let's remove the c-spine translational DoFs from the experimental 
        % motion tracking terms
        q_i = Xkj(1:2:end, i + 1);
        
        
        if ~isempty(special_tracking_names)
            q_special_diff = q_i(special_indices_ind, :) - Qs_track_kj(special_indices_ind, i + 1); 
            q_special_term = q_special_term + W.q_special * B(i + 1) * J_special(q_special_diff) * mesh_T;

        end
        
        q_new_indices = setdiff(local_tracking_indices, special_indices_ind)';
        q_diff = q_i(q_new_indices, :) - Qs_track_kj(indices_to_track, i + 1);
        q_term = q_term + W.q * B(i + 1) * J_q(q_diff) * mesh_T;

        % term for special DoFs
        if ~isempty(no_tracking_names)
            q_no_tracking_term = q_no_tracking_term + W.q_no_tracking * B(i + 1) * J_no_tracking(q_i(no_tracking_indices_ind', :)) * mesh_T;

        end

        % Q dot
        q_dot_i = Xkj(2:2:end, i + 1);
        q_dot_diff = q_dot_i(local_tracking_indices, :) - Qdots_track_kj(:, i + 1);
        q_dot_term = q_dot_term + W.q_dot * B(i + 1) * J_q_dot(q_dot_diff) * mesh_T;


        % Q dot dot: Accelerations
        acc_term = acc_term + W.acc * B(i + 1) * J_acc(Aj(:, i)) * mesh_T;


        % GRF
        if ~isempty(grf_path)
            Ti_GRF_scaled = Ti(num_q_all + dll_grf_indices, 1)./scaling_GRF';
            GRF_diff = Ti_GRF_scaled - GRF_track_j(:, i);
            GRF_term = GRF_term + W.GRF * B(i + 1) * J_GRF(GRF_diff) * mesh_T;
        end


        % residual terms
        for n = 1:length(res_keys)
            key   = res_keys{n};
            J_func = J_res.(key);
            res_term.(key) = res_term.(key) + W.residuals * B(i + 1) * J_func(res_sym.(key)(:, i)) * mesh_T;

            residuals_term = residuals_term + res_term.(key);

        end

 
               
        % ------------ add them up ----------- %
        J = q_term + q_special_term + q_dot_term + acc_term + q_no_tracking_term + ...
            GRF_term + residuals_term;
        
        
        % --------------------------------------------------------------- %
        %                      Equality constraints                       %
        % --------------------------------------------------------------- %        
        
        % Rigid body dynamics: 
        % 
        % here, we want to impose the constraint
        % that the derivative of the positional data at the current
        % collocation point corresponds to the q_dot at the same
        % collocation point. To do this, we use the C matrix, 
        % from the collocation scheme, to compute the
        % derivative approximation at the collocation points of the
        % segment.This is because, technically, IPOPT treats the
        % two discretised trajectories (q and q_dot)
        % individually.Therefore, we need to impose this equality
        % constraint to make sure that the physics is correct.
        
        Q_nsc_dot  = Xkj_nsc(1:2:end, :) * C(:, i + 1);  % [num_q, d + 1] @ [d+1, 1] -> [num_q, 1]
        Qdots_nsc_dot  = Xkj_nsc(2:2:end, :) * C(:, i + 1);  % [num_q, 1] @ [d + 1, 1] -> [num_q, 1]    
        Qdotj_nsc = Xkj_nsc(2:2:end, i + 1); % velocity
        eq_constr{end+1} = (mesh_T * Qdotj_nsc - Q_nsc_dot)./x_scaling.Qs(:, q_ind_indices)';
        eq_constr{end+1} = (mesh_T * Aj_nsc(:, i) - Qdots_nsc_dot)./x_scaling.Qsdot(:, q_ind_indices)';
      

        % torque activation dynamics: 
        % 
        % we impose the constraint that the derivative of the torque
        % activation at the current collocation point must be equal to the
        % analytical derivative that is computed via the activation
        % dynamics formula.
        a_a_dot  = a_akj * C(:, i + 1); % [num_actuators, d + 1] @ [d+1, 1] -> [num_actuators, 1]

        % torque activation dynamics (explicit formulation)   
        da_dt_i = activation_dynamics_function(e_ak, a_akj(:, i + 1));
        eq_constr{end+1} = (mesh_T * da_dt_i - a_a_dot);
        

 
        % Path constraints 
        % --------------------------------------------------------------- %
        % here, we want to impose the constraint that 
        % Computed torque from CasADi should be equal to the net moments
        % coming out of the OpenSim model. We do this only for the DoFs
        % that aren't spanned by the muscles.
        Ti_torques = Ti(actuator_indices, 1)./actuator_sf;
        eq_constr{end+1} = Ti_torques - a_akj(:, i + 1);
        
        % ----- residuals ----- %
        for n = 1:length(res_keys)
            key   = res_keys{n};
            res_tau = Ti(res_sf_indices.(key), 1)./res_sf_value.(key);
            eq_constr{end+1} = res_sym.(key)(:, i) - res_tau;

        end


    end
    
    
    fprintf('Collocation loop finished\n');



    eq_constr = vertcat(eq_constr{:});

    % Now we define a CasADi function that takes in the design variables at
    % the collocation points and outputs the cost function (J) and the sets
    % of constraints.    
    func_out =  {eq_constr, J, q_term, q_dot_term, acc_term};

    f_coll = Function('f_coll', func_in, func_out);

    % save expression graph
    f_coll.save(char(fullfile(output_dir, 'f_J.casadi')));
    
    
    % register function as parallel form across the number of segments of
    % the trajectory.
    f_coll_map = f_coll.map(N, parallelMode, num_threads);


    % finally, we evaluate everything that was built symbolically
    [eq_constr_all, J_all, q_term_all,...
        q_dot_term_all, acc_term_all] = f_coll_map(func_map_in{:});

    fprintf('Function has been evaluated\n');
    


    % ------------------------------------------------------------------- %
    % add constraints to opti struct
    opti.subject_to(eq_constr_all == 0);

    lbg = opti.value(opti.lbg);
    ubg = opti.value(opti.ubg);
    n_eq   = sum(lbg == ubg);
    n_ineq = sum(lbg ~= ubg);
    
    fprintf('Equality constraints before continuity : %d\n', n_eq);
    fprintf('Inequality constraints: %d\n', n_ineq);
    fprintf('Decision variables    : %d\n', opti.nx);
    
    % --------------------------------------------------------------- %
    %                       mesh end points                           %
    % --------------------------------------------------------------- %
    % Loop over segments
    Q_mesh = X(1:2:end, :);
    Q_col = X_col(1:2:end, :);
    
    Qdot_mesh = X(2:2:end, :);
    Qdot_col = X_col(2:2:end, :);
    
    col_indices = [1, 2, 3];
    for k=1:N

        q_kj = [Q_mesh(:, k), Q_col(:, col_indices)];
        qdot_kj = [Qdot_mesh(:, k), Qdot_col(:, col_indices)];
        a_akj = [a_a(:,k), a_a_col(:, col_indices)];

        % Add equality constraints (next interval starts with end values of 
        % states from previous interval)
        opti.subject_to(Q_mesh(:, k + 1) == q_kj * D);
        opti.subject_to(Qdot_mesh(:, k + 1) == qdot_kj * D);
        opti.subject_to(a_a(:, k + 1) == a_akj * D);

        % update collocation indices
        col_indices = col_indices + 3;
    end

    lbg = opti.value(opti.lbg);
    ubg = opti.value(opti.ubg);
    n_eq   = sum(lbg == ubg);

    fprintf('Equality constraints after continuity : %d\n', n_eq);
    
    % sum objective function over all mesh segments
    J_sum = sum(J_all);

    % --------------------------------------------------------------- %
    %                          NLP solver                             %
    % --------------------------------------------------------------- %           
    opti.minimize(J_sum);
    options.ipopt.hessian_approximation = 'limited-memory';
    options.ipopt.mu_strategy  = 'adaptive';
    options.ipopt.max_iter = config_struct.("optimiser").("max_iters");
    tolerance = config_struct.("optimiser").("tolerance");
    options.ipopt.tol = 1*10^(-tolerance);
    options.ipopt.print_timing_statistics = 'yes';
    options.ipopt.nlp_scaling_method = 'none';
    options.ipopt.obj_scaling_factor = 1;
    % options.ipopt.linear_solver = 'mumps';          % <-- add this
    % options.ipopt.mumps_mem_percent = 200;          % <-- and this
    opti.solver('ipopt', options);

    % --------------------------------------------------------------- %
    %                          Solve problem                          %
    % --------------------------------------------------------------- %   
    [w_opt, stats, g_opt, lambda_x, lambda_g] = solve_NLPSOL(opti, options);  


end
