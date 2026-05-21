function res = track_sim(trial_id, trial_dir, dll_filepath)

% --------------------------------------------------------------------------
% track_sim
%   This function aims at tracking experimental data while solving an
%   Optimal Control Problem.

% INPUTs:
%    - trial_id (str): the ID of the specific trial under investigation

%    - trial_dir (str): path to the specific trial directory

%   - dll_filepath (str): path to the dll specific to the current trial 
%     being tracked.


% The function assumes the following structure in trial_dir:
% /trial_dir
%  |__.osim file  --> scaled OpenSim model
%  |__/grf   --> where the GRF .mot files are contained
%  |__/ik    --> where the IK .mot files are contained


%
% --------------------------------------------------------------------------
    import casadi.*
    import org.opensim.modeling.*
    
    % read project config file
    json_str = fileread("config.json");
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

    pathMuscleModel = [pwd,'/muscle_model'];
    addpath(genpath(pathMuscleModel));   

    pathCasADiFunctions = [pwd,'/casadi_functions'];
    addpath(genpath(pathCasADiFunctions));

    pathBounds = [pwd,'/bounds'];
    addpath(genpath(pathBounds));

    pathGetters = [pwd,'/getters'];
    addpath(genpath(pathGetters));

    pathPolynomial = [pwd,'/polynomials'];
    addpath(genpath(pathPolynomial));
    
    
    % Collocation scheme
    N = config_struct.("collocation").("number_of_segments");   % number of mesh intervals
    d = config_struct.("collocation").("num_points"); % number of collocation points per mesh interval
    method = config_struct.("collocation").("method"); % collocation method
    
    [tau_root, C, D, B] = collocationScheme(d, method); % collocation scheme.
    
    
    % Load external functions
    F = external('F', dll_filepath); 

    % -------------------------- OpenSim Model -------------------------- % 
    model_filepath = fullfile(trial_dir, "/scaled_model.osim");
    model = Model(model_filepath);

     % ---------- sets ----------- %
    coordinateSet = model.getCoordinateSet();
    forceSet = model.getForceSet();
    muscleSet = forceSet.getMuscles();
    
    
    % ---------- info ----------- %
    q_names = getItemNames(coordinateSet);
    muscle_names = getItemNames(muscleSet);
    MTparameters = getMTparameters(model, muscle_names);

    % generic parameters
    num_body_dof = 6;  % the number of theoretical DoFs that a rigid body has
    num_q_all = coordinateSet.getSize();  % number of coordinates (independent and dependent)
    num_muscles = muscleSet.getSize();  % number of muscles



    % ------------------------ Experimental Data ------------------------ % 
    
    % ------ IK ------ %
    ik_filepath = fullfile(trial_dir, "/ik/", trial_id + "_ik.mot");
    Qs = readMotQs(ik_filepath, 20);

    dof_indices_all = 1:num_q_all;

    % read from config structure the names of the DoFs that are spanned by
    % the muscles.
    dof_names = config_struct.("dof_names_spanned_by_muscles");
    
    % find which indices of Qs correspond to the DoFs that are spanned by
    % the neck muscles.
    dof_indices = cellfun(@(name) find(strcmp(q_names, name)), dof_names);

    % now find the indices of all the other DoFs that aren't spanned by
    % muscles.
    other_indices = setdiff(dof_indices_all, dof_indices);

    % we register the number of the DoFs that aren't spanned by the muscles
    % as "number of actuators". This simply means that the force actuator
    % for these coordinates is modelled as an idealised actuator, and it
    % doesn't have all the inherent complexity of a muscle.
    num_actuators = size(other_indices, 1);

    % read names of dependent coordinates
    dependent_coord_names = config_struct.("dependent_coord_names");
    % find dependent coordinates indices
    dependent_coord_idx = cellfun(@(name) find(strcmp(q_names, name)), dependent_coord_names);
    
    
    % differentiate between number of independent and dependent coords 
    num_q_dep = size(dependent_coord_names, 1);
    num_q_ind = num_q_all - num_q_dep;
    
    
    % load Ground Reaction Forces
    grf_filepath = fullfile(trial_dir, "/grf/", trial_id + ".mot");
    GRFs = readMotGrf(grf_filepath, 20);
    experimental_force_indices = config_struct.("experimental_force_indices");
    grf_indices = [experimental_force_indices.rGRF'; experimental_force_indices.lGRF'];
    grm_indices = [experimental_force_indices.rGRM'; experimental_force_indices.lGRM'];


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
    Qs.allinterpfilt = interp1(Qs.time, Qs.allfilt, time_intervals);
    
    % find Qs values at each collocation point along the trajectory
    Qs.allinterpfilt_col = interp1(Qs.time, Qs.allfilt, time_grid');


    % --- GRF --- %
    % find indices where the 2 items of time_opt are 
    hz_GRF = 2000;  % make this more flexible, probably computing it from the grf.mot file as in the IK case is ideal
    dt_GRF = 1 / hz_GRF;
    time_expi.GRF(1) = find((GRF.time<(time_opt(1) + dt_GRF/2)) & (GRF.time>=(time_opt(1) - dt_GRF/2)));
    time_expi.GRF(2) = find((GRF.time<(time_opt(2) + dt_GRF/2)) & (GRF.time>=(time_opt(2) - dt_GRF/2)));

    % ----------------------------- Bounds  ----------------------------- %
    [bounds, scaling] = getBounds(Qs, GRFs, num_q, num_muscles, num_act, grf_indices, grm_indices);

    
    % -------------------------- Initial Guess  ------------------------- %
    guess = getGuess(Qs, num_q, num_muscles, num_act, scaling);
    
    
    % ------------------- Experimental Data Scaling  ------------------- %
    % we have the scaled experimental data stored in the 'guess' struct.
   
    Qs_scaled = guess.Qs_all(:, 1:2:end);
    Qs_scaled_col = guess.Qs_col(:, 2:2:end);
    
    Qdots_scaled = guess.Qs_all(:, 1:2:end);
    Qdots_scaled_col = guess.Qs_col(:, 2:2:end);
    
    
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
    opti.subject_to(bounds.X.lower' < X_col < bounds.X.upper');
    % initial condition
    opti.set_initial(X_col, guess.Qs_col');
    
    % 2) mesh end-points
    points = N + 1;  
    X = opti.variable(dims, points);
    opti.subject_to(bounds.X.lower' < X < bounds.X.upper');
    opti.set_initial(X, guess.Qs_all');

    % ---------- Muscles ---------- %
    % Activation
    dims = num_muscles;
    points = d * N;
    
    % 1) collocation points
    a_col = opti.variable(dims, points);
    opti.subject_to(bounds.a.lower'< a_col < bounds.a.upper');
    opti.set_initial(a_col, guess.a_col');

    % 2) mesh end-points
    points = N;
    a = opti.variable(dims, points);
    opti.subject_to(bounds.a.lower'< a < bounds.a.upper');
    opti.set_initial(a, guess.a'); 

    % MT force
    dims = num_muscles;
    points = d * N;
    
    % 1) collocation points
    FTtilde_col = opti.variable(dims, points);
    opti.subject_to(bounds.FTtilde.lower'< FTtilde_col < bounds.FTtilde.upper');
    opti.set_initial(FTtilde_col, guess.FTtilde_col');

    % 2) mesh end-points
    points = N;
    FTtilde = opti.variable(dims, points);
    opti.subject_to(bounds.FTtilde.lower'< FTtilde < bounds.FTtilde.upper');
    opti.set_initial(FTtilde, guess.FTtilde');


    
    % ----- Torque Actuators ----- %
    dims = num_act;
    points = d * N;
    
    % 1) collocation points
    a_a_col = opti.variable(dims, points);
    opti.subject_to(bounds.a_a.lower'< a_a_col < bounds.a_a.upper');
    opti.set_initial(a_a_col, guess.a_a_col');

    % 2) mesh end-points
    points = N;
    a_a = opti.variable(dims, points);
    opti.subject_to(bounds.a_a.lower'< a_a < bounds.a_a.upper');
    opti.set_initial(a_a, guess.a_a');  


    % ----------------------- Controls  ----------------------- %
    % ----- Muscles ----- %
    % Time derivative of muscle activations (states) at mesh points
    dims = num_muscles;
    points = N;
    
    % end points
    vA = opti.variable(num_muscles, N);
    opti.subject_to(bounds.vA.lower'*ones(1,N) < vA < bounds.vA.upper'*ones(1,N));
    opti.set_initial(vA, guess.vA'); 

    % ----- Actuator Excitation ----- %
    e_a = opti.variable(nq.arms, N);
    opti.subject_to(bounds.e_a.lower'*ones(1, N) < e_a < bounds.e_a.upper'*ones(1,N));
    opti.set_initial(e_a, guess.e_a');

    % Define "slack" controls
    % Time derivative of muscle-tendon forces (states) at collocation points
    dFTtilde_col = opti.variable(NMuscle, d * N);
    opti.subject_to(bounds.dFTtilde.lower'*ones(1,d*N) < dFTtilde_col < ...
            bounds.dFTtilde.upper'*ones(1, d * N));
    opti.set_initial(dFTtilde_col, guess.dFTtilde_col');
    
    % Time derivative of Qdots (states) at collocation points
    A_col = opti.variable(nq.all, d * N);
    opti.subject_to(bounds.Qdotdots.lower'*ones(1, d * N) < A_col < ...
            bounds.Qdotdots.upper'*ones(1, d * N));
    opti.set_initial(A_col, guess.Qdotdots_col'); 


    % ------------------------------------------------------------------- %
    % The following section uses the "MX" object from the CasADi library
    % to register symbolic variables that will later be used in CasADi 
    % "Function" objects, which dictate how these variables interact with
    % each other.

    % The suffix "k" indicates the variable at the mesh end-points.
    % "j", on the other hand, at the collocation points.
    % ------------------------------------------------------------------- %

    % ---------- state-vector ---------- %
    % activation
    ak = MX.sym('ak', num_muscles);
    aj = MX.sym('akmesh', num_muscles, d);
    akj = [ak aj];
    
    % muscle-tendon unit force
    FTtildek = MX.sym('FTtildek', num_muscles); 
    FTtildej = MX.sym('FTtildej', num_muscles, d);
    FTtildekj = [FTtildek FTtildej];
    
    % q and q_dot
    Xk = MX.sym('Xk', 2 * num_q);   % shape: [num_q * 2]
    Xj = MX.sym('Xj', 2 * num_q, d); % shape: [num_q * 2, d]
    Xkj = [Xk Xj];  % shape: [num_q * 2, d + 1]
    
    a_ak = MX.sym('a_ak', num_act);
    a_aj = MX.sym('a_akmesh', num_act, d);
    a_akj = [a_ak a_aj];
    
    
    % ---------- controls ---------- %
    % da/dt
    vAk = MX.sym('vAk', num_muscles);
    
    % torque actuator excitation
    e_ak = MX.sym('e_ak', num_act);
    
    % dF/dt
    dFTtildej = MX.sym('dFTtildej', NMuscle,d);
    
    % q_dotdot: remember, we're using implicit formulation. Therefore,
    % accelerations are treated as "controls".
    Aj = MX.sym('Aj', num_q, d); 

    % ----------Pelvis residuals ---------- %
    pelvis_res_j = MX.sym('pelvis_res_j', num_body_dof, d);


    % ------- Experimental Data to Track ------- %
    Qs_scaled_k = MX.sym('Qs_scaled_k', size(Qs_scaled, 2));
    Qs_scaled_j = MX.sym('Qs_scaled_j', size(Qs_scaled, 2), d);
    Qs_scaled_kj = [Qs_scaled_k Qs_scaled_j];
    
    Qdots_scaled_k = MX.sym('Qdots_scaled_k', size(Qs_scaled, 2));
    Qdots_scaled_j = MX.sym('Qdots_scaled_j', size(Qs_scaled, 2), d);
    Qdots_scaled_kj = [Qdots_scaled_k Qdots_scaled_j];
    
    num_grfs = 6;
    GRF_scaled_k = MX.sym('GRF_scaled_k', num_grfs); 
    GRF_scaled_j = MX.sym('GRF_scaled_j', num_grfs, d); 
    GRF_scaled_kj = [GRF_scaled_k GRF_scaled_j];

    num_grm = 6;
    GRM_scaled_k = MX.sym('GRM_scaled_k', num_grm); 
    GRM_scaled_j = MX.sym('GRM_scaled_j', num_grm, d); 
    GRM_scaled_kj = [GRM_scaled_k GRM_scaled_j];


    % initialise set of constraint vector
    eq_constr = {}; % equality constraint vector
    ineq_constr1 = {}; % inequality constraint
    ineq_constr2 = {}; 
    ineq_constr3 = {}; 
    ineq_constr4 = {}; 
    g_names_coll = {}; % Initialize names of constraints at collocation points
    
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
    Xkj_nsc = Xkj.*(scaling.QsQdots');  % shape: [num_q * 2, d + 1]
     
    FTtildekj_nsc = FTtildekj.*(scaling.FTtilde');
    
    dFTtildej_nsc = dFTtildej.*scaling.dFTtilde;
    Aj_nsc = Aj.*(scaling.Qdotdots');  
    
    vAk_nsc = vAk.*scaling.vA;  

     
    % ---------- CasADi functions ---------- %
    f_muscle = buildMuscleFunction(pathMuscleModel);
    force_equilibrium = buildForceEquilibriumFunc(pathMuscleModel, num_muscles, MTparameters);
    
    % Torque actuation dynamics
    activation_dynamics_function = torque_activation_dynamics_casadi(pathMuscleModel, num_act);

    % --- sum of squares --- %
    % q
    J_q = sum_of_squares('q', num_q);
    
    % q_dot
    J_q_dot = sum_of_squares('q_dot', num_q);

    % q_dot_dot
    J_q_dot_dot = sum_of_squares('q_dot_dot', num_q);

    % muscle activation
    J_muscles_act = sum_of_squares('muscle_act', num_muscles);

    % da/dt
    J_muscles_act_der = sum_of_squares('muscle_act_der', num_muscles);

    % MT unit force
    J_MT_unit = sum_of_squares('MT_unit', num_muscles);

    % GRFs
    J_GRF = sum_of_squares('GRF', num_grfs);

    % GRMs
    J_GRM = sum_of_squares('GRM', num_grm);

    % Pelvis residuals
    J_pelvis = sum_of_squares('pelvis', num_body_dof);

    
    % --- joint moments --- %
    
    % active elements
    muscle_spanning_joint_file = [polyResultsPath "/muscle_spanning_joint_INFO_subject_GC.mat"];
    muscle_spanning_joint_info = load(muscle_spanning_joint_file);
    
    % loop through the DoFs spanned by the model muscles 
    for n = 1:size(dof_indices, 1)
        % retrieve row of 0's and 1's for the current DoF
        dof_values = muscle_spanning_joint_info(:, n);
        
        % find which muscles span the current DoF
        indices = find(dof_values);

        % get number of muscles that span the current DoF
        D = size(indices, 1);

        % get DoF name
        dof_name = dof_names(:, n);

        % generate function
        M_functions.(dof_name) = compute_active_moment(D);
        
    end

    % passive elements
    M_passive = compute_passive_moment();

    
    % Initialise cost function value: J
    J = 0;

    % Read the weights of the cost function from config file
    W = config_struct.("W");
    
    % for one segment:
    % loop through the collocation points
    for i = 1:d 
        
        [Xkj_nsc_all, Aj_all] = apply_constraints(Xkj_nsc(:, i + 1), Aj(:, i + 1));
        
        % evaluate external function
        Ti = F([Xkj_nsc_all, Aj_all]);

        % --- compute: lMT (muscle-tendon length), d_lMT/dt, MA(moment arm)
        
        % isolate DoFs that are spanned by muscles
        Q_kj_nsc = Xkj_nsc(1:2:end, i + 1); % +1 because first point is the mesh beginning
        Qdot_kj_nsc = Xkj_nsc(2:2:end, i + 1);
        
        dof_i = Q_kj_nsc(dof_indices);  
        dof_dot_i = Qdot_kj_nsc(dof_indices);

        [lMT_i, vMT_i, dM_i] = f_muscle(dof_i, dof_dot_i);

        % Get muscle-tendon forces and derive Hill-equilibrium 
        [hill_err_i, FT_i, ~, ~, ~] = force_equilibrium(akj(:, i+1), ...
            FTtildekj_nsc(:, i+1), dFTtildej_nsc(:, i), lMT_i, vMT_i);

        
        % --------------------------------------------------------------- %
        % Use the C matrix, from the collocation scheme, to compute the
        % derivative approximation at the collocation points of the
        % segment.
        % --------------------------------------------------------------- %
        Q_nsc_dot  = Xkj_nsc(1:2:end,:) * C(:, i+1);  % [num_q, d + 1] @ [d+1, d] -> [num_q, d]
        Qdots_nsc_dot  = Xkj_nsc(2:2:end,:) * C(:, i+1);  % [num_q, d + 1] @ [d+1, d] -> [num_q, d]          
        
        FTtilde_nsc_dot  = FTtildekj_nsc * C(:, i+1);% [num_q, d + 1] @ [d+1, d] -> [num_q, d]
   
        a_dot  = akj * C(:, i+1); % [num_muscle, d + 1] @ [d+1, d] -> [num_muscles, d]
        
        a_a_dot  = a_akj * C(:, i+1); % [num_actuators, d + 1] @ [d+1, d] -> [num_actuators, d]
        
        % ---------------- Cost Function ---------------- %
        % difference between optimised and observed Qs
        q_diff = Q_kj_nsc - Qs_scaled_kj(:, i + 1);
        q_dot_diff = Qdots_nsc_dot - Qdot_kj_nsc(:, i + 1);
        

        % ---------------- Terms ---------------- %
        % Q
        q_term = W.q * B(i + 1) * J_q(q_diff) * mesh_T;

        % Q dot
        q_dot_term = W.q_dot * B(i + 1) * J_q_dot(q_dot_diff) * mesh_T;

        % GRF
        Ti_GRF_scaled = Tj(num_q + 1, 1)./scaling.GRF';
        GRF_term = W.GRF * B(i + 1) * (J_GRF(Ti_GRF_scaled - GRF_scaled_kj(:, i))) * mesh_T;

        % GRM
        Ti_GRM_scaled = Tj(GRFi.all, 1)./scaling.GRM';
        GRM_term = W.GRM * B(i + 1) * (J_GRM(Ti_GRM_scaled - GRM_scaled_kj(:, i))) * mesh_T;

        % muscle activation
        act_term = W.a * B(i + 1) * (J_muscles_act(akj(:, i + 1))) * mesh_T;

        % muscle activation time derivative
        act_der_term = W.vA * B(i + 1) * (J_muscles_act_der(vAk) * mesh_T);

        % Q dot dot: Accelerations
        q_dot_dot_term = W.acc * B(i + 1) * (J_q_dot_dot(Aj(:, i)) * mesh_T);

        % derivative of MT force
        dMTf_dt_term = B(i + 1) * (J_MT_unit(dFTtildej(:, i)) * mesh_T);

        % pelvis residual term
        pelvis_term = W.pelvis * B(i + 1) * (J_pelvis(pelvis_res_j(:, i)) * mesh_T);
        
        
        % ------------ add them up ----------- %
        J = q_term + q_dot_term + q_dot_dot_term + GRF_term + GRM_term + ...
            act_term + act_der_term + dMTf_dt_term + pelvis_term;
        
        
        % --------------------------------------------------------------- %
        %                      Equality constraints                       %
        % --------------------------------------------------------------- %        
        
        % Muscle activation time derivative
        eq_constr{end+1} = (mesh_T * vAk_nsc - a_dot)./scaling.a;
        
        % Contraction dynamics (implicit formulation)     
        eq_constr{end+1} = (mesh_T * dFTtildej_nsc(:, i) - FTtilde_nsc_dot)./scaling.FTtilde';
        
        % Skeleton dynamics (implicit formulation)               
        Qdotj_nsc = Xkj_nsc(2:2:end, i +1 ); % velocity
        eq_constr{end+1} = (mesh_T * Qdotj_nsc - Q_nsc_dot)./scaling.QsQdots(1:2:end)';
        eq_constr{end+1} = (mesh_T * Aj_nsc(:, i) - Qdots_nsc_dot)./scaling.QsQdots(2:2:end)';


        % Torque actuator activation dynamics (explicit formulation)   
        da_dt_i = activation_dynamics_function(e_ak, a_akj(:, i+1));
        eq_constr{end+1} = (mesh_T * da_dt_i - a_a_dot)./scaling.a_a;

        % Computed torque from CasADi should be equal to the net moments
        % coming out of the OpenSim model. We do this only for the DoFs
        % that aren't spanned by the muscles.
        eq_constr{end+1} = Ti(other_indices, 1)./scaling.a_a - a_akj(:, i+1);

        % Path constraints 
        % --------------------------------------------------------------- %
        % here, we want to impose the constraint that the net joint moments 
        % coming out of the external function must be equal to the sum of:
        % the moment cause by muscles and passive joint elements.

        % loop through the DoFs spanned by the model muscles 
        for n = 1:size(dof_indices, 1)
            % retrieve row of 0's and 1's for the current DoF
            dof_values = muscle_spanning_joint_info(:, n);
        
            % find which muscles span the current DoF
            indices = dof_values;

            % get DoF name
            dof_name = dof_names(:, n);

            % find the index of the DoF wrt the full list of coordinates
            idx = strcmp(q_names, dependent_coord_name);


            % ----------- active elements ----------- %
            % get moment arms
            moment_arms = dM_i(indices, n);

            % get forces
            forces = FT_i(indices);

            % retrieve function
            M_function = M_functions.(dof_name);
            M_computed = M_function(moment_arms, forces);
            
            % add difference to equality constraint
            eq_constr{end+1} = Ti(idx, 1) - M_computed;            
        
        end

        % Contraction dynamics (implicit formulation)
        eq_constr{end+1} = Hilldiffj;

        
        % --------------------------------------------------------------- %
        %                      Inequality constraints                     %
        % --------------------------------------------------------------- %
        % Activation dynamics (implicit formulation)   
        tact = 0.015;
        tdeact = 0.06;
        act1 = vAk_nsc + akj(:, i + 1)./(ones(size(akj(:, i + 1),1),1) * tdeact);
        act2 = vAk_nsc + ak./(ones(size(ak,1),1)*tact);
        ineq_constr1{end+1} = act1;
        ineq_constr2{end+1} = act2; 
        

    end

    eq_constr = vertcat(eq_constr{:});
    ineq_constr1 = vertcat(ineq_constr1{:});
    ineq_constr2 = vertcat(ineq_constr2{:});
    ineq_constr3 = vertcat(ineq_constr3{:});

    % Now we define a CasADi function that takes in the design variables at
    % the collocation points and outputs the cost function (J) and the sets
    % of constraints.
    f_coll = Function('f_coll', {Xk, Xj, Aj, ak, aj, vAk, ...
        FTtildek, FTtildej, dFTtildej, ...
        a_ak, a_aj, e_ak, dFTtildej,...
        Qs_scaled_k, Qs_scaled_j, ...
        Qdots_scaled_k,Qdots_scaled_j,...
        GRF_scaled_k, GRM_scaled_j, ...
        GRM_scaled_k, GRM_scaled_j},...
        {eq_constr, ineq_constr1, ineq_constr2, ineq_constr3, J});

    % register function as parallel form across the number of segments of
    % the trajectory.
    f_coll_map = f_coll.map(N, parallelMode, num_threads);


    % finally, we evaluate everything that was built symbolically
    [eq_constr_num, ineq_constr1_num, ineq_constr2_num, ...
        coll_ineq_constr3, J_num] = f_coll_map(X(:,1:end-1), X_col, A_col, ...
        a(:, 1:end-1), a_col, vA,...
        FTtilde(:, 1:end-1), FTtilde_col, dFTtilde_col, ...
        a_a, a_a_col,e_a,...
        ...
        Qs_scaled(1:end-1,:)', ...
        Qs_scaled_col',...
        Qdots_scaled(1:end-1,:)',...
        Qdots_scaled_col',...
        GRF.val.allinterp_col(:, 2:end)',...
        GRF.MorGF.allinterp_col(:, 2:end)');

end
