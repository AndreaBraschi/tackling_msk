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
    
    solveProblem = true;            % Set true to solve the optimal control problem.
    analyseResults = true;          % Set true to analyze the results.
    loadResults = false;            % Set true to load the results of the optimization.
    saveResults = true;             % Set true to save the results of the optimization.
    checkBoundsIG = false;          % visualize guess-bounds
    writeMotionFiles = true;        % Set true to write motion files for use in OpenSim GUI
    saveOptimalTrajectories = true; % Set true to save optimal trajectories
    writeIKmotion=true;             % Set true to write .mot file


    % parallelisation settings
    parallelMode = 'thread';
    NThreads = 8; % Number of threads used in parallel.

    
    % --- add paths of subdirectories where we'll be picking functions from --- % 
    pathmain = pwd;
    [pathRepo,~,~] = fileparts(pathmain);
    
    pathSettings = [pathRepo,'/Settings'];
    addpath(genpath(pathSettings));

    % various utility functions
    pathUtils = [pathRepo,'/utils'];
    addpath(genpath(pathUtils));

    % collocation scheme
    pathCollocationScheme = [pathRepo,'/collocationScheme'];
    addpath(genpath(pathCollocationScheme));

    pathMuscleModel = [pathRepo,'/MuscleModel'];
    addpath(genpath(pathMuscleModel));   

    pathExternalFunctions = [pathRepo,'/externalFunctions/'];
    if exist(pathExternalFunctions, 'dir')
       cd(pathExternalFunctions);
    else
        warning('Cannot change directory: %s does not exist.', pathExternalFunctions);

    end
    
    pathBounds = [pathRepo,'/Bounds'];
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
    [tau_root, C, D, B] = collocationScheme(d, method); % collocation scheme. See function to understand how state variables and their derivatives are computed at the collocation points
    
    
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

    % 4) let's get some information on the model
    num_q = coordinateSet.getSize();  % number of coordinates
    num_muscles = muscleSet.getSize();  % number of muscles
    num_act = actuatorSet.getSize();  % number of actuators

    
% -------------------------- Experimental Data -------------------------- % 
    % load IK
    nametrial.GRF   = [nametrial.id, "grf"];
    nametrial.IK    = [nametrial.id, "_IK"];

    Qs = readMotQs(fullfile(trial_path, nametrial.IK));
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
    
    
    
    % ------------------- Experimental Data Scaling  ------------------- %
    %% TODO: 
    % this is a 3rd repetition of the same interpolation process. It
    % happens during 'getBounds' and 'getGuesses'. We could simply compute
    % the spline coefficients once from the main script and pass that as
    % input to the other functions where interpolation is performed.
    % 
    Qs_scaled = (Qs.allinterpfilt)/scaling.Qs;
    Qs_scaled_col = interp1(interval(1:N+1), Qs_res_interpfilt_scaled_aux, time_grid);
    
    for i=2:size(Qs.allinterpfilt,2)
        Q_spline(i-1)=spline(Qs.allinterpfilt(:,1),Qs.allinterpfilt(:,i));
        Qdot_spline(i-1)=fnder(Q_spline(i-1),1);
        Qdots.allinterpfilt(:,1)=Qs.allinterpfilt(:,1);
        Qdots.allinterpfilt(:,i)=ppval(Qdot_spline(i-1),Qdots.allinterpfilt(:,1));
    end
    Qdots_res_interpfilt_scaled=Qdots.allinterpfilt(:,2:end)./scaling.Qdots;
    Qdots_res_interpfilt_scaled_aux=Qdots_res_interpfilt_scaled(:,[1:4 6:14 20:34]); %exclude knee internal dofs
    Qdots_res_interpfilt_scaled_aux_col=interp1(interval(1:N+1),Qdots_res_interpfilt_scaled_aux,tgrid_col);
    
    
    

    % -------------------------- Initial Guess  ------------------------- %
    guess = getGuess(Qs, num_q, num_muscles, num_act, scaling);
    
    

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
    points = d*N;
    
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
    opti.subject_to(bounds.e_a.lower'*ones(1,N) < e_a < bounds.e_a.upper'*ones(1,N));
    opti.set_initial(e_a, guess.e_a');

    % Define "slack" controls
    % Time derivative of muscle-tendon forces (states) at collocation points
    dFTtilde_col = opti.variable(NMuscle, d*N);
    opti.subject_to(bounds.dFTtilde.lower'*ones(1,d*N) < dFTtilde_col < ...
            bounds.dFTtilde.upper'*ones(1,d*N));
    opti.set_initial(dFTtilde_col, guess.dFTtilde_col');
    
    % Time derivative of Qdots (states) at collocation points
    A_col = opti.variable(nq.all, d*N);
    opti.subject_to(bounds.Qdotdots.lower'*ones(1,d*N) < A_col < ...
            bounds.Qdotdots.upper'*ones(1,d*N));
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
    Xk = MX.sym('Xk', 2 * num_q);
    Xj = MX.sym('Xj', 2 * num_q, d);
    Xkj = [Xk Xj];
    
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
    Aj = MX.sym('Aj', nq.all,d); 



    % ------- Experimental Data to Track ------- %
    Qs_scaled_k = MX.sym('Qs_scaled_k', size(Qs_res_interpfilt_scaled_aux, 2));
    Qs_scaled_j = MX.sym('Qs_scaled_j',size(Qs_res_interpfilt_scaled_aux, 2), d);
    Qs_scaled_kj = [Qs_scaled_k Qs_scaled_j];
    
    Qdots_scaled_k = MX.sym('Qdots_scaled_k', size(Qs_res_interpfilt_scaled_aux, 2));
    Qdots_scaled_j = MX.sym('Qdots_scaled_j', size(Qs_res_interpfilt_scaled_aux, 2), d);
    Qdots_scaled_kj = [Qdots_scaled_k Qdots_scaled_j];
    
    num_grfs = 6;
    GRF_scaled_k = MX.sym('GRF_scaled_k', num_grfs); 
    GRF_scaled_j = MX.sym('GRF_scaled_j', num_grfs, d); 
    GRF_scaled_kj = [GRF_scaled_k GRF_scaled_j];

    num_grm = 6;
    GRM_scaled_k = MX.sym('GRM_scaled_k', num_grm); 
    GRM_scaled_j = MX.sym('GRM_scaled_j', num_grm, d); 
    GRM_scaled_kj = [GRM_scaled_k GRM_scaled_j];



end