function res = track_sim(model, trial_path, dll_filename, useReducedPolynomials, err_poly, Options, W)

% --------------------------------------------------------------------------
% track_sim
%   This function aims at tracking experimental data while solving an
%   Optimal Control Problem.

% INPUTs:
%   - model: scaled OpenSim model file of the form Model(file.osim).
% 
%   - trial_id (str): id of the trial that will be processed

%   - dll_filename (str): dll file that contains the external function used
%     in the optimisation.
% 
%   - useReducedPolynomials
%   
%   - err_poly:
%
%   - Options:
%
%   - W
%
% OUTPUT:
%   - guess -
%   * initial guess values for all optimisation variables
% 
% --------------------------------------------------------------------------


    import casadi.*
    import org.opensim.modeling.*
    
    % parallelisation settings
    parallelMode = 'thread';
    num_threads = 8; % Number of threads used in parallel.

    
    % --- add paths of subdirectories where we'll be picking functions from --- % 
    pathmain = pwd;
    [pathRepo, ~, ~] = fileparts(pathmain);
    
    pathSettings = [pathRepo,'/Settings'];
    addpath(genpath(pathSettings));

    % various utility functions
    pathUtils = [pathRepo,'/utils'];
    addpath(genpath(pathUtils));

    % collocation scheme
    pathCollocationScheme = [pathRepo,'/collocationScheme'];
    addpath(genpath(pathCollocationScheme));

    pathMuscleModel = [pathRepo,'/muscle_model'];
    addpath(genpath(pathMuscleModel));   

    pathCasADiFunctions = [pathRepo,'/casadi_functions'];
    addpath(genpath(pathSettings));


    pathExternalFunctions = [pathRepo,'/externalFunctions/'];
    if exist(pathExternalFunctions, 'dir')
       cd(pathExternalFunctions);
    else
        warning('Cannot change directory: %s does not exist.', pathExternalFunctions);

    end
    
    pathBounds = [pathRepo,'/bounds'];
    addpath(genpath(pathBounds));

    pathPolynomial = [pathRepo,'/Polynomials_GC'];
    addpath(genpath(pathPolynomial));
    
    
   
    nametrial.id = trial_id;

    % simulation settings - add them to Options structure
    Options.usePelvisResMom = 1;    % use pelvis residual moments

    Options.useReducedPolynomials = useReducedPolynomials;  % use reduced polynomial coeffcients
    Options.useReoptimizedPoly = 1;

    Options.err_poly = err_poly;
    Options.maxsmoothness = 'MellowMax'; %options: logSum, MellowMax, nosmooth

    tol_ipopt = 4;    % tolerance (means 1e-(tol_ipopt))

    setup.derivatives = 'AD_Recorder'; % Algorithmic differentiation / Recorder     

    % Available linear solvers
    linear_solvers = {'mumps','ma27','ma57','ma77','ma86','ma97'}; 
    if Options.useReducedPolynomials
        if Options.useReoptimizedPoly
            poly=['red' num2str(Options.err_poly) '_reopt'];
        else
            poly=['red' num2str(Options.err_poly)];
        end
    else
        poly=['full' num2str(Options.err_poly)];
    end
    
    % The filename used to save the results depends on the settings 
    if ~exist('savename_suffix', 'var')
        savename_suffix = '';
    end
    savename = ['_', nametrial.id, '_', num2str(Options.IGn), '_poly', poly, savename_suffix];
    savename2 = ['_', nametrial.id, '_', num2str(Options.IGn), '_poly', poly];
    
    
    % Collocation scheme
    N = 40;   % number of mesh intervals
    d = 3; % number of collocation points per mesh interval
    method = 'radau'; % collocation method
    [tau_root, C, D, B] = collocationScheme(d, method); % collocation scheme.
    
    
    % Load external functions
    % The external function performs inverse dynamics through the
    % OpenSim/Simbody C++ API. This external function is compiled as a dll from
    % which we create a Function instance using CasADi in MATLAB. 
    % We use different external functions. A first external function extracts 
    % several parameters of the bodies to which the contact spheres are attached.
    % The contact forces are then computed in MATLAB and are inputs of the
    % second external function in which the skeleton dynamics is described. The
    % motivation for this decoupling is to limit the number of times we need to
    % build the model. By defining the contact model in MATLAB, we only need to
    % build the model once per external function, whereas keeping the contact
    % model in the external function would require re-building the model during
    % the optimization.
    dll_path = [pathExternalFunctions, dll_filename]; 
    F = external('F', dll_path); 
    cd(pathmain);


    % --------------------- Model Information --------------------- % 
    % 1) Mass
    state = model.initSystem();
    bodyMass = model.getTotalMass(state);
    bodyWeight = bodyMass * 9.81;

    % 2) Model sets
    coordinateSet = model.getCoordinateSet();
    forceSet = model.getForceSet();
    muscleSet = forceSet.getMuscles();
    actuatorSet = forceSet.getActuators();
    q_names = getItemNames(coordinateSet);


    % 3) Muscle Tendon Unit parameters
    muscleNames = getItemNames(muscleSet);
    MTparameters = getMTparameters(model, muscleNames);

    % 4) information on the model
    num_body_dof = 6;  % the number of theoretical DoFs that a rigid body has
    num_q = coordinateSet.getSize();  % number of coordinates
    num_muscles = muscleSet.getSize();  % number of muscles
    num_act = actuatorSet.getSize();  % number of actuators

    
% -------------------------- Experimental Data -------------------------- % 
    % load IK
    nametrial.GRF   = [nametrial.id, "grf"];
    nametrial.IK    = [nametrial.id, "_IK"];

    Qs = readMotQs(fullfile(trial_path, nametrial.IK));
    
    % find which indices of Qs correspond to the DoFs that are spanned by
    % the neck muscles.
    dof_names = {'arm_flex_l', 'arm_add_l', 'arm_rot_l', 'arm_flex_r', 'arm_add_r', 'arm_rot_r', ...
    'auxt1jnt_r3', 'auxt1jnt_r1', 'auxt1jnt_r2', ...
    'auxt1jnt_t3', 'auxt1jnt_t1', 'auxt1jnt_t2', ...
    'aux7jnt_t3', 'aux7jnt_t1', 'aux7jnt_t2', ...
    'aux6jnt_t3', 'aux6jnt_t1', 'aux6jnt_t2', ...
    'aux5jnt_t3', 'aux5jnt_t1', 'aux5jnt_t2', ...
    'aux4jnt_t3', 'aux4jnt_t1', 'aux4jnt_t2', ...
    'aux3jnt_t3', 'aux3jnt_t1', 'aux3jnt_t2', ...
    'aux2jnt_r3', 'aux2jnt_r1', 'aux2jnt_r2'};
    
    dof_indices = cellfun(@(name) find(strcmp(Qs.colheaders, name)), dof_names, 'UniformOutput', false);

    dt_ik = Qs.time(2) - Qs.time(1);

    % load Ground Reaction Forces
    GRFs = readMotGrf(fullfile(trial_path, nametrial.GRF));


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
    Qs.allinterpfilt = interp1(Qs.allfilt(:, 1), Qs.allfilt, time_intervals);
    
    % find Qs values at each collocation point along the trajectory
    Qs.allinterpfilt_col = interp1(Qs.allfilt(:,1), Qs.allfilt, time_grid');


    % --- GRF --- %
    % find indices where the 2 items of time_opt are 
    hz_GRF = 2000;  % make this more flexible, probably computing it from the grf.mot file as in the IK case is ideal
    dt_GRF = 1 / hz_GRF;
    time_expi.GRF(1) = find((GRF.time<(time_opt(1) + dt_GRF/2)) & (GRF.time>=(time_opt(1) - dt_GRF/2)));
    time_expi.GRF(2) = find((GRF.time<(time_opt(2) + dt_GRF/2)) & (GRF.time>=(time_opt(2) - dt_GRF/2)));

    % ----------------------------- Bounds  ----------------------------- %
    [bounds, scaling] = getBounds(Qs, GRFs, num_q, num_muscles, num_act);

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
    dims = 2 * num_q;
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
    
    % Structure W with 
    W.Qs = 125;
    
    
    % ---------- CasADi functions ---------- %
    f_muscle = buildMuscleFunction(pathMuscleModel);
    force_equilibrium = buildForceEquilibriumFunc(pathMuscleModel, num_muscles, MTparameters);
    
    % Torque actuation dynamics
    activation_dynamics_function = torque_activation_dynamics_casadi(pathMuscleModel, num_act);

    % --- sum of squares --- %
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
    
    % for one segment:
    % loop through the collocation points
    for i = 1:d 
        [Xkj_nsc_all, Aj_all] = apply_constraints(Xkj_nsc(:, i + 1), Aj(:, i + 1));
        Ti = F([Xkj_nsc_all, Aj_all]);

        % --- compute: lMT (muscle-tendon length), d_lMT/dt, MA(moment arm)
        
        % isolate DoFs that are spanned by muscles
        Q_kj_nsc = Xkj_nsc(1:2:end, i + 1); % +1 because first point is the mesh beginning
        Qdot_kj_nsc = Xkj_nsc(2:2:end, i + 1);
        
        dof_i = Q_kj_nsc(dof_indices);  
        dof_dot_i = Qdot_kj_nsc(dof_indices);

        [lMT_i, vMT_i, dM_i] = f_muscle(dof_i, dof_dot_i);

        % Get muscle-tendon forces and derive Hill-equilibrium 
        [hill_err_i, FT_i, ~, ~, ~] = force_equilibrium(akj(:, i+1), FTtildekj_nsc(:, i+1), dFTtildej_nsc(:, i), lMT_i, vMT_i);

        
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
        
        % vectorise weight of the error
        weight_qs_vec = W.Qs * ones(num_q, 1);

        % ---------------- Terms ---------------- %
        % Q
        q_term = B(i+1) * (q_diff' * diag(weight_qs_vec) * q_diff) * mesh_T;

        % Q dot
        q_dot_term = B(i + 1) * J_q_dot(q_dot_diff);

        % GRF
        Ti_GRF_scaled = Tj(GRFi.all, 1)./scaling.GRF';
        GRF_term = B(i + 1) * (J_GRF(Ti_GRF_scaled - GRF_scaled_kj(:, i))) * mesh_T;

        % GRM
        Ti_GRM_scaled = Tj(GRFi.all, 1)./scaling.GRM';
        GRM_term = B(i + 1) * (J_GRM(Ti_GRM_scaled - GRM_scaled_kj(:, i))) * mesh_T;

        % muscle activation
        act_term = B(i + 1) * (J_muscles_act(akj(:, i + 1))) * mesh_T;

        % muscle activation time derivative
        act_der_term = B(i + 1) * (J_muscles_act_der(vAk) * mesh_T);

        % Q dot dot: Accelerations
        q_dot_dot_term = B(i + 1) * (J_q_dot_dot(Aj(:, i)) * mesh_T);

        % derivative of MT force
        dMTf_dt_term = B(i + 1) * (J_MT_unit(dFTtildej(:, i)) * mesh_T);

        % pelvis residual term
        pelvis_term = B(i + 1) * (J_pelvis(pelvis_res_j(:, i)) * mesh_T);
        
        
        % ------------ add them up ----------- %
        J = q_term + q_dot_term + q_dot_dot_term + GRF_term + GRM_term + ...
            act_term + act_der_term + dMTf_dt_term + pelvis_term;
        
        
        % --------------------------------------------------------------- %
        %   compute (scaled) equality constraints at collocation points   %
        % --------------------------------------------------------------- %
        
        eq_constr{end+1} = (mesh_T * vAk_nsc - a_dot)./scaling.a;
        
        % Contraction dynamics (implicit formulation)     
        eq_constr{end+1} = (mesh_T * dFTtildej_nsc(:, i) - FTtilde_nsc_dot)./scaling.FTtilde';
        
        % Skeleton dynamics (implicit formulation)               
        Qdotj_nsc = Xkj_nsc(2:2:end, i +1 ); % velocity
        eq_constr{end+1} = (mesh_T * Qdotj_nsc - Q_nsc_dot)./scaling.QsQdots(1:2:end)';
        eq_constr{end+1} = (mesh_T * Aj_nsc(:, i) - Qdots_nsc_dot)./scaling.QsQdots(2:2:end)';


        % Arm activation dynamics (explicit formulation)   
        da_dt_i = activation_dynamics_function(e_ak, a_akj(:, i+1));
        eq_constr{end+1} = (mesh_T * da_dt_i - a_a_dot)./scaling.a_a;


        % --------------------------------------------------------------- %
        %                       Path constraints                          %
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

            % get moment arms
            moment_arms = dM_i(indices, n);

            % get forces
            forces = FT_i(indices);

            % retrieve function
            M_function = M_functions.(dof_name);
            M_computed = M_function(moment_arms, forces);

            % add difference to equality constraint
            eq_constr{end+1} = Ti(jointi.knee_add.r,1)-(T_knee_add_r + Tau_passj.knee_add.r);            
        
        end
        
    end



        


 
    end



end