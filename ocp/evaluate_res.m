function evaluate_res(trial_id, trial_dir, dll_filepath)

    load checkpoint.mat w_opt stats;
    load scaling.mat scaling;
    load guess.mat guess;


    import casadi.*
    import org.opensim.modeling.*
    
    % read project config file
    json_str = fileread("config.json");
    config_struct = jsondecode(json_str);
    
    % --- add paths of subdirectories where we'll be picking functions from --- %     
    % various utility functions
    pathUtils = [pwd,'/utils'];
    addpath(genpath(pathUtils));

    % collocation scheme
    pathCollocationScheme = [pwd,'/collocationScheme'];
    addpath(genpath(pathCollocationScheme));

    pathMuscleModel = [pwd,'/muscle_model'];
    addpath(genpath(pathMuscleModel));   

    pathMusclePoly = [pwd,'/muscle_polynomials'];
    addpath(genpath(pathMusclePoly));   


    pathCasADiFunctions = [pwd,'/numerical_functions'];
    addpath(genpath(pathCasADiFunctions));

    pathBounds = [pwd,'/bounds'];
    addpath(genpath(pathBounds));

    pathGetters = [pwd,'/getters'];
    addpath(genpath(pathGetters));

    pathPolynomial = [pwd,'/muscle_polynomials'];
    addpath(genpath(pathPolynomial));

    pathConstraints = [pwd,'/kinematic_coupling'];
    addpath(genpath(pathConstraints));


   
    % Collocation scheme
    N = config_struct.("collocation").("number_of_segments");   % number of mesh levals
    d = config_struct.("collocation").("num_points"); % number of collocation points per mesh interval
    method = config_struct.("collocation").("method"); % collocation method
    
    [tau_root, C, D, B] = collocationScheme(d, method); % collocation scheme.

    
    % Load external functions
    F = external('F', dll_filepath); 
    % load information relevant to external function from config file
    dll_force_indices = config_struct.("dll_force_indices");
    dll_grf_indices = [dll_force_indices.rGRF', dll_force_indices.lGRF'];
    dll_grm_indices = [dll_force_indices.rGRM', dll_force_indices.lGRM'];

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
    % read experimental qs indices values
    pelvis_y_idx = config_struct.("experimental_q_indices").("pelvis_y");

    % import scaled factors
    scale_factors = config_struct.("scale_factors");
    pelvis_y = scale_factors.("pelvis_y");

    
    % ------ IK ------ %
    ik_filepath = fullfile(trial_dir, "/ik/", trial_id + "_ik.mot");
    Qs = readMotQs(ik_filepath, 15);
    Qs.allfilt(:, pelvis_y_idx + 1) = Qs.allfilt(:, pelvis_y_idx + 1) + ...
        Qs.allfilt(:, pelvis_y_idx + 1) * pelvis_y;

    % find the indices of the c-spine translational DoFs
    spine_t_names = config_struct.("c_spine_translational_dof_names");
    spine_t_indices = cellfun(@(name) find(strcmp(q_names, name)), spine_t_names);

    % zero the values of Qs corresponding to the spine translational DoFs.
    % We can do it here as we're not going to track the experimental data,
    % and we're going to make our life easier when it comes to guess
    % generation (0) and bounds ([-1 1])
    Qs.allfilt(:, spine_t_indices + 1) = 0;

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


    % read names of dependent coordinates
    q_dep_names = config_struct.("dependent_coord_names");
    % find dependent coordinates indices
    q_dep_indices = cellfun(@(name) find(strcmp(q_names, name)), q_dep_names);


    % independent coordinates indices
    q_ind_indices = setdiff(dof_indices_all, q_dep_indices);
    q_ind_names = q_names(q_ind_indices);

    % differentiate between number of independent and dependent coords 
    num_q_dep = size(q_dep_names, 1);
    num_q_ind = num_q_all - num_q_dep;  

    % now, let's find the indices of the c-spine translational DoFs within
    % the set of independent DoFs
    spine_t_ind_indices = cellfun(@(name) find(strcmp(q_ind_names, name)), spine_t_names');
    local_ind_indices = 1:num_q_ind;
    local_tracking_indices = setdiff(local_ind_indices, spine_t_ind_indices);
    
    
    % now let's get the indices of the q's to track
    q_tracking_indices = setdiff(q_ind_indices, spine_t_indices');
    num_q_tracking = size(q_tracking_indices, 2);

    independent_act_indices = setdiff(other_indices, q_dep_indices');
    independent_act_indices = [independent_act_indices(7:end-6), independent_act_indices(end-2:end)];
    num_actuators = size(independent_act_indices, 2);

    % load Ground Reaction Forces
    grf_filepath = fullfile(trial_dir, "/grf/", trial_id + ".mot");
    GRFs = readMotGrf(grf_filepath, 20);
    experimental_force_indices = config_struct.("experimental_force_indices");
    grf_indices = [experimental_force_indices.rGRF', experimental_force_indices.lGRF'];
    grm_indices = [experimental_force_indices.rGRM', experimental_force_indices.lGRM'];

    
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


    % q
    Qs_scaled_all = guess.Qs_all(:, 1:2:end);
    Qs_scaled_col_all = guess.Qs_col(:, 1:2:end);
        
    % retrieve just independent coordinates
    Qs_scaled = Qs_scaled_all(:, q_tracking_indices);
    Qs_scaled_col = Qs_scaled_col_all(:, q_tracking_indices)';
    
    % q dot
    Qdots_scaled_all = guess.Qs_all(:, 2:2:end);
    Qdots_scaled_col_all = guess.Qs_col(:, 2:2:end);
    % retrieve just independent coordinates
    Qdots_scaled = Qdots_scaled_all(:, q_tracking_indices);
    Qdots_scaled_col = Qdots_scaled_col_all(:, q_tracking_indices)';
    

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
    GRM_col = interp1(GRFs.time, GRFs.data(:, grm_indices), time_grid');
    
    % GRF;
    GRF_scaled = (GRF_col./scaling.GRF)';


    % GRM;
    GRM_scaled = (GRM_col./scaling.GRM)';

    

    % --------------- Retrieve q's -------------------- %
    % collocation points
    dims = 2 * num_q_ind;
    points = N * d;  % we have d-collocation points x N-segments
    X_col = zeros(dims, points);
    u = 1;
    for n = 1:points
        X_col(:, n) = w_opt(u:u + dims - 1, :);
        u = u + dims;
    end

    % mesh end-points
    points = N + 1;
    X = zeros(dims, points);
    for n = 1:points
       X(:, n) = w_opt(u:u + dims - 1, :);
       u = u + dims;
    end


    % --------------- Retrieve muscle activation -------------------- %
    dims = num_muscles;
    
    % collocation points
    points = d * N;
    a_col = zeros(dims, points);
    for n = 1:points
       a_col(:, n) = w_opt(u:u + dims - 1, :);
       u = u + dims;
    end


    % 2) mesh end-points
    points = N + 1;
    a = zeros(dims, points);
    for n = 1:points
       a(:, n) = w_opt(u:u + dims - 1, :);
       u = u + dims;
    end

    % --------------- Retrieve MT force -------------------- %
    dims = num_muscles;
    
    % 1) collocation points
    points = d * N;
    FTtilde_col = zeros(dims, points);
    for n = 1:points
       FTtilde_col(:, n) = w_opt(u:u + dims - 1, :);
       u = u + dims;
    end

    % 2) mesh end-points
    points = N + 1;
    FTtilde = zeros(dims, points);
    for n = 1:points
       FTtilde(:, n) = w_opt(u:u + dims - 1, :);
       u = u + dims;
    end


    % --------------------- Retrieve torque actuators ---------------------- %
    dims = num_actuators;
    points = d * N;
    
    % 1) collocation points
    a_a_col = zeros(dims, points);
    for n = 1:points
       a_a_col(:, n) = w_opt(u:u + dims - 1, :);
       u = u + dims;
    end
    
    
    % 2) mesh end-points
    points = N + 1;
    a_a = zeros(dims, points);
    for n = 1:points
       a_a(:, n) = w_opt(u:u + dims - 1, :);
       u = u + dims;
    end



    % --------------------- Retrieve muscle vA ------------------------ %
    dims = num_muscles;
    points = N + 1;
    
    % end points
    vA = zeros(dims, points);
    for n = 1:points
       vA(:, n) = w_opt(u:u + dims - 1, :);
       u = u + dims;
    end

    
    % ------------------ Retrieve actuator excitation ------------------- %
    dims = num_actuators;
    e_a = zeros(dims, points);
    for n = 1:points
       e_a(:, n) = w_opt(u:u + dims - 1, :);
       u = u + dims;
    end



    % ---------------- Retrieve MT force time derivative ---------------- %
    dims = num_muscles;
    points = d * N;
    dFTtilde_col = zeros(dims, points);
    for n = 1:points
       dFTtilde_col(:, n) = w_opt(u:u + dims - 1, :);
       u = u + dims;
    end
    
    
    % ------------------- Retrieve accelerations ------------------- %
    dims = num_q_ind;
    points = d * N;
    A_col = zeros(dims, points);
    for n = 1:points
       A_col(:, n) = w_opt(u:u + dims - 1, :);
       u = u + dims;
    end
    
    
    
    % ----------------------- Retrieve residuals  ----------------------- %
    % pelvis
    dims = num_body_dof;
    points = d * N;
    pelvis_res_col = zeros(dims, d * N);
    for n = 1:points
       pelvis_res_col(:, n) = w_opt(u:u + dims - 1, :);
       u = u + dims;
    end
    
    
    % punching bag
    dims = 3;
    pb_res_col = zeros(dims, d * N);
    for n = 1:points
       pb_res_col(:, n) = w_opt(u:u + dims - 1, :);
       u = u + dims;
    end


    % we need to make use of the unscaled version of the data
    X_nsc = zeros(size(X_col));
    X_nsc(1:2:end, :) = X_col(1:2:end,:).*scaling.Qs(:, q_ind_indices)';
    X_nsc(2:2:end, :) = X_col(2:2:end,:).*scaling.Qsdot(:, q_ind_indices)';


    % plot kinematics
    interleaved = reshape(repmat(q_ind_names, 2, 1), 1, []);
    for index = 1:size(interleaved, 2)

        figure('Visible', 'off');
        plot(time_grid, X_col(index, :), 'b-', 'LineWidth', 1.5, 'DisplayName', 'IPOPT');
        hold on;
        plot(time_grid, guess.Qs_col_ind(:, index), 'r--', 'LineWidth', 1.5, 'DisplayName', 'Experimental');
        hold off;

        xlabel('Time (s)');
        ylabel('Coordinate value');

        if mod(index, 2) ~= 0
            label = interleaved{index};
        else
            label = [interleaved{index} 'Vel'];
        end

        title(label);


        legend('show');
        grid on;

        saveas(gcf, sprintf('C:/Users/ab3758/Documents/PhD/msk/ocp/figures/qs/adjusted_pelvis_ty/coord_%02d_%s.png', index, interleaved{index}));
        close(gcf);
    end

    A_nsc = A_col.*(scaling.Qsdotdot(:, q_ind_indices)');  
   
    % initialise GRFs
    num_grfs = 6;
    GRFs = zeros(num_grfs, d * N);
    GRMs = zeros(num_grfs, d * N);
        
    % Read the weights of the cost function from config file
    W = config_struct.("W");

    
    % loop through the collocation points
    col_points = d * N;
    
    % initialise cost function terms
    q_term = zeros(col_points, 1);
    q_dot_term = zeros(col_points, 1);
    GRF_term = zeros(col_points, 1);
    GRM_term = zeros(col_points, 1);
    act_term = zeros(col_points, 1);
    act_der_term = zeros(col_points, 1);
    q_dot_dot_term = zeros(col_points, 1);
    dMTf_dt_term = zeros(col_points, 1);
    pelvis_term = zeros(col_points, 1);
    pb_term = zeros(col_points, 1);
    
    
    u = 1;
    for i = 1:col_points 
        
        % current q-part of state vector
        x_i = X_nsc(:, i);

        % current q accelerations
        acc_i = A_nsc(:, i);

        [x_all_i, acc_all_i] = apply_constraints(x_i, acc_i, q_names, q_ind_indices);

        % evaluate external function
        Ti = F(vertcat(x_all_i, acc_all_i));
        GRFs(:, i) = full(Ti(num_q_all + dll_grf_indices, 1));
        GRMs(:, i) = full(Ti(num_q_all + dll_grm_indices, 1));


        % ---------------- Cost Function Terms ---------------- %
        % Q
        % let's remove the c-spine translational DoFs from the experimental 
        % motion tracking terms
        q_i = X_col(1:2:end, i);
        q_diff = q_i(local_tracking_indices, :) - Qs_scaled_col(:, i);
        q_term(i, :) = W.q * B(u + 1) * sum(q_diff.^2) * mesh_T;

        % Q dot
        q_dot_i = X_col(2:2:end, i);
        q_dot_diff = q_dot_i(local_tracking_indices, :) - Qdots_scaled_col(:, i);
        q_dot_term(i, :) = W.q_dot * B(u + 1) * sum(q_dot_diff.^2) * mesh_T;

        % GRF
        Ti_GRF_scaled = GRFs(:, i)./scaling.GRF';
        GRF_diff = Ti_GRF_scaled - GRF_scaled(:, i);
        GRF_term(i, :) = W.GRF * B(u + 1) * sum(GRF_diff.^2) * mesh_T;

        % GRM
        Ti_GRM_scaled = GRMs(:, i)./scaling.GRM';
        GRM_diff = Ti_GRM_scaled - GRM_scaled(:, i);
        GRM_term(i, :) = W.GRM * B(u + 1) * sum(GRM_diff.^2) * mesh_T;

        % muscle activation
        act_term(i, :) = W.a * B(u + 1) * (sum(a_col(:, i).^2)) * mesh_T;

        % muscle activation time derivative
        act_der_term(i, :) = W.vA * B(u + 1) * sum(vA(:).^2) * mesh_T;

        % Q dot dot: Accelerations
        q_dot_dot_term(i, :) = W.acc * B(u + 1) * sum(A_col(:, i).^2) * mesh_T;

        % derivative of MT force
        dMTf_dt_term(i, :) = W.u * B(u + 1) * sum(dFTtilde_col(:, i).^2) * mesh_T;

        % pelvis residual term
        pelvis_term(i, :) = W.pelvis * B(u + 1) * sum(pelvis_res_col(:, i).^2) * mesh_T;

        % punching bag residual term
        pb_term(i, :) = W.pelvis * B(u + 1) * sum(pb_res_col(:, i).^2) * mesh_T;

        if mod(i, 3) == 0
           u = 1;
        else
           u = u + 1;
        end
        
           
    end
    

    J = sum(q_term) + sum(q_dot_term) + sum(GRF_term) + sum(GRM_term) + ...
        sum(act_term) + sum(act_der_term) + sum(q_dot_dot_term) + ...
        sum(dMTf_dt_term) + sum(pelvis_term) + sum(pb_term);

    % sum each term across collocation points
    term_values = [
        sum(q_term), ...
        sum(q_dot_term), ...
        sum(GRF_term), ...
        sum(GRM_term), ...
        sum(act_term), ...
        sum(act_der_term), ...
        sum(q_dot_dot_term), ...
        sum(dMTf_dt_term), ...
        sum(pelvis_term), ...
        sum(pb_term)
    ];

    term_labels = {
        'q tracking', ...
        'q dot tracking', ...
        'GRF', ...
        'GRM', ...
        'muscle activation', ...
        'activation derivative', ...
        'q dot dot', ...
        'dMTf/dt', ...
        'pelvis', ...
        'pb'
    };

    figure;
    bar(term_values);
    set(gca, 'XTick', 1:length(term_labels), ...
             'XTickLabel', term_labels, ...
             'XTickLabelRotation', 45);
    ylabel('Cost Contribution');
    title('Cost Function Component Contributions');
    grid on;
    saveas(gcf, sprintf('C:/Users/ab3758/Documents/PhD/msk/ocp/figures/qs/adjusted_pelvis_ty/cost_function_histogram.png'));
    close(gcf);
        

end

