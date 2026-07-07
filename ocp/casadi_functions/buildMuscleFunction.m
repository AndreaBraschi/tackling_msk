function f_muscle = buildMuscleFunction(pathMusclePoly)

import casadi.*
muscle_spanning_joint_file = fullfile(pathMusclePoly, "muscle_spanning_joints_info_modified.mat");
muscle_spanning_joint_info = load(muscle_spanning_joint_file);

MuscleInfo_file = fullfile(pathMusclePoly, "muscle_info_opt.mat");
MuscleInfo = load(MuscleInfo_file);

muscle_joint_info_fieldnames = fieldnames(muscle_spanning_joint_info);
muscle_joint_info_field = muscle_joint_info_fieldnames{1};

num_q = size(muscle_spanning_joint_info.(muscle_joint_info_field), 2);
num_muscles = size(muscle_spanning_joint_info.(muscle_joint_info_field), 1);

qin     = SX.sym('qin', 1, num_q);
qdotin  = SX.sym('qdotin', 1, num_q);
lMT     = SX(num_muscles, 1);
vMT     = SX(num_muscles, 1);
dM      = SX(num_muscles, num_q);


% loop through muscles
for i=1:num_muscles     
    
    % find which DoFs the current muscle spans
    dof_indices = find(muscle_spanning_joint_info.(muscle_joint_info_field)(i, :) == 1); 
    
    % retrieve polynomial order: i.e. 3rd
    order = MuscleInfo.MuscleInfo.muscle(i).order;

    % compute muscle length and moment arm, symbolically, given the set of
    % q's that the current muscle spans.
    [mat, diff_mat_q] = mvpoly_sym(qin(1, dof_indices), order);

    % find the indices where the coefficients are non-zero
    idx_nonzeros = find(MuscleInfo.MuscleInfo.muscle(i).coeff{order});
    
    lMT(i, 1) = mat(idx_nonzeros) * MuscleInfo.MuscleInfo.muscle(i).coeff{order}(idx_nonzeros);
    vMT(i, 1) = 0;
    dM(i, 1:num_q) = 0;
    nr_dof_crossing = length(dof_indices); 
    
    for dof_nr = 1:nr_dof_crossing

        % moment arm
        dM(i, dof_indices(dof_nr)) = (-(diff_mat_q(idx_nonzeros, dof_nr)))' * MuscleInfo.MuscleInfo.muscle(i).coeff{order}(idx_nonzeros);
        
        % rate of change of MT length: d_lMT/dt
        vMT(i,1) = vMT(i,1) + (-dM(i, dof_indices(dof_nr))* qdotin(1, dof_indices(dof_nr)));
            
     end 
 end

 f_muscle = Function('f_muscle', {qin, qdotin}, {lMT,vMT,dM});

end