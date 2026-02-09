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

    
    % --- add paths of the repo where we'll be picking functions from --- % 
    pathmain = pwd;
    [pathRepo,~,~] = fileparts(pathmain);
    
    pathSettings = [pathRepo,'/Settings'];
    addpath(genpath(pathSettings));

    pathFileFunctions = [pathRepo,'/filesFunctions'];
    addpath(genpath(pathFileFunctions));

    pathCollocationScheme = [pathRepo,'/CollocationScheme'];
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
    Options.maxsmoothness='MellowMax'; %options: logSum, MellowMax, nosmooth

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
    [tau_root, C, D, B] = CollocationScheme(d, method); % collocation scheme. See function to understand how state variables and their derivatives are computed at the collocation points
    
    
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
    dll_path = [pathExternalFunctions, dll_name]; 
    F = external('F', dll_path); 
    cd(pathmain);




    % --------------------- Model Information --------------------- % 
    % 1) Mass
    state = model.initSystem();
    bodyMass = model.getTotalMass(state);
    bodyWeight = bodyMass * 9.81;

    % 2) Coordinates
    coordinateSet = model.getCoordinateSet();
    q_names = getItemNames(coordinateSet);


    % 3) Muscle Tendon Unit parameters
    forceSet = model.getForceSet();
    muscleSet = forceSet.getMuscles();
    muscleNames = getItemNames(muscleSet);
    MTparameters = getMTparameters(model, muscleNames);

    % 4) muscles actuating the independed DoFs

    
% -------------------------- Experimental Data -------------------------- % 
    % load IK
    nametrial.GRF   = [nametrial.id, "grf"];
    nametrial.IK    = [nametrial.id, "_IK"];

    Qs = readMotFile(fullfile(trial_path, nametrial.IK));
    dt_ik = Qs.time(2) - Qs.time(1);

    % load Ground Reaction Forces
    GRFs = readMotGrf(fullfile(trial_path, nametrial.GRF));


    % read initial and final time from IK 
    time_opt = [Qs.time(1, 1) Qs.time(end, 1)];


    % ------------ Interpolation ------------ %
    % compute the time window of each mesh
    step = (time_opt(2) - time_opt(1)) / N;

    % create a time vector that goes from time_opt(1) to time_opt(2)
    % .window by window.
    intervals = time_opt(1):step:time_opt(2);
    % now, let's create a vector that contains time stamps that correspond
    % to every single collocation point along the trajectory.
    time_grid = zeros(N*d,1);
    time_grid(1:d:end) = intervals(1:end-1) + tau_root(2)*step; % 1st collocation point of every mesh
    time_grid(2:d:end) = intervals(1:end-1) + tau_root(3)*step; % 2nd collocation point of every mesh
    time_grid(3:d:end) = intervals(1:end-1) + tau_root(4)*step; % 3rd collocation point of every mesh


    % --------- IK --------- %

    % --------- GRF --------- %
    % find indices where the 2 items of time_opt are 
    hz_GRF = 2000;  % make this more flexible, probably computing it from the grf.mot file as in the IK case is ideal
    dt_GRF = 1 / hz_GRF;
    time_expi.GRF(1) = find((GRF.time<(time_opt(1) + dt_GRF/2)) & (GRF.time>=(time_opt(1) - dt_GRF/2)));
    time_expi.GRF(2) = find((GRF.time<(time_opt(2) + dt_GRF/2)) & (GRF.time>=(time_opt(2) - dt_GRF/2)));


    % CasADi function --> what does this do? 



end