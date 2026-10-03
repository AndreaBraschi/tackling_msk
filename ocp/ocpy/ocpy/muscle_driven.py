import json
import numpy as np
from scipy.interpolate import CubicSpline
import casadi as ca
import opensim as osim
from typing import Optional, List
from osim_utils.read import readStoFile

# ----------- local imports ---------- #
from .collocation.schemes import radau
from .getters import *
from .utils.filters import butterworth

from .bounds import *
from .guess import *


from .casadi_functions import torque_activation_dynamics_casadi
from .casadi_functions import sum_of_squares
from .casadi_functions.apply_constraints import apply_constraints

from .opt import solve_NLP

def muscle_driven(model_path: str, ik_path: str, dll_path: str, config_filepath: str, kinematic_coupling_path: str,
                      output_dir: str, grf_path: Optional[str] = None):


    # ---------- read config file ----------- #
    with open(config_filepath) as f:
        cfg = json.load(f)


    #------------ load external function (rigid body dynamics) ----------- #
    F = ca.external('F', dll_path)
    dll_force_idx = cfg["dll_force_indices"]
    dll_grf_idx = np.array(dll_force_idx["rGRF"] + dll_force_idx["lGRF"]) - 1  # 0-based

    # -------------------------- OpenSim Model -------------------------- #
    model = osim.Model(model_path)

    # ---------- sets ----------- #
    coordinateSet = model.getCoordinateSet()

    # ---------- model info ---------- #
    # get coordinate names - all: independent and dependent
    q_names = get_item_names(coordinateSet)
    num_q_all = len(q_names)
    q_all_indices = range(num_q_all)

    # read names of dependent coordinates
    q_dep_names = cfg["dependent_coord_names"]

    # find dependent coordinates indices
    q_dep_indices = [i for i, q_name in enumerate(q_names) if q_name in q_dep_names]
    q_dep_names_check = np.array(q_names)[q_dep_indices]

    # find independednt coordinate indices and names
    q_ind_indices = np.array([i for i, q_name in enumerate(q_names) if q_name not in q_dep_names])
    q_ind_names = np.array(q_names)[q_ind_indices]

    # differentiate between number of independent and dependent coords
    num_q_dep = len(q_dep_names)
    num_q_ind = num_q_all - num_q_dep

    # check if the user doesn't want specific coordinates not to be tracked
    no_tracking_names = cfg["tracking"]["no_tracking"]

    if no_tracking_names:
        # find indices of listed names
        no_tracking_indices = [i for i, q_name in enumerate(q_names) if q_name in no_tracking_names]
        # how many they are
        num_no_tracking_indices = len(no_tracking_indices)

        # find global indices of the coordinates that are to be tracked
        q_tracking_indices = [i for i in q_ind_indices if i not in no_tracking_indices]

        # now, this can be a little confusing at this point of the code. However, normally we only track the independent
        # coordinates from experimental data. The independent coordinate set represents a subset of indices. In this
        # new subset which goes from 0:end-1, we need to find the "local" indices of the names that the user provided.
        no_tracking_indices_l = [i for i, name in enumerate(q_ind_names) if name in no_tracking_names]
        tracking_indices_l = [i for i, name in enumerate(q_ind_names) if name not in no_tracking_names]

    else:
        no_tracking_indices = []
        local_ind_indices = range(num_q_ind)
        tracking_indices_l = local_ind_indices
        q_tracking_indices = q_ind_indices

    num_q_tracking = len(q_tracking_indices)
    tracking_indices = range(num_q_tracking)

    #  check if the user has declared a specific q to be tracked differently (i.e. different weight: higher or lower)
    special_tracking_names = cfg["tracking"]["special_tracking"]
    if special_tracking_names:
        special_q_indices = [i for i, q_name in enumerate(q_names) if q_name in special_tracking_names]
        num_special_indices = len(special_q_indices)

        # now, we do the same thing as before: we must express the indices also in the "local" subset.
        # let's find the indices of the special DoFs within the set of independent DoFs
        special_q_indices_l = [i for i, name in enumerate(q_ind_names) if name in special_tracking_names]

        # then we need a list of indices that will allow us to extract the coordinates that will be assigned a global
        # tracking weight
        indices_to_track = [i for i in tracking_indices if i not in special_q_indices_l]

    else:
        indices_to_track = tracking_indices


    # ------------------------ Experimental Data ------------------------ #
    # ---------- IK ---------- #
    Qs_raw = readStoFile(ik_path)
    time_raw = Qs_raw["time"].values
    dt_raw = time_raw[1] - time_raw[0]
    ikCoordinates = Qs_raw.values[:, 1:].copy()   # the .copy() command makes sure that ikCoordinates isn't just a view of the DataFrame and can be modified in place

    # filter
    order = 2
    cutoff = 20

    # convert to radians
    # 1) check which indices contain translational DoFs
    trans_indices = np.array([i for i, col in enumerate(Qs_raw.columns) if any(col.endswith(s) for s in ['_tx', '_ty', '_tz',
                                                                                                '_t3', '_t1', '_t2'])]) - 1

    trans_mask = np.zeros(ikCoordinates.shape[-1], dtype=bool)
    trans_mask[trans_indices] = True

    # 2) convert degrees --> radians
    ikCoordinates[:, ~trans_mask] *= (np.pi / 180)

    Qs = butterworth(time_raw, ikCoordinates.T, dt_raw, order, cutoff, 'low')

    #  assign the initial value prescribed by the user
    if no_tracking_names:
       Qs[:, no_tracking_indices] = 0

    # see if the user provided any GRF
    if grf_path:
        GRFs_raw = readStoFile(grf_path)
        GRFs = butterworth(GRFs_raw["time"].values, GRFs_raw.values[:, 1:].T, dt_raw, order, cutoff, 'low')
        dt_GRF = GRFs_raw["time"][1] - GRFs_raw["time"][0]
        grf_indices = np.array(cfg["experimental_force_indices"]["rGRF"] + cfg["experimental_force_indices"]["lGRF"]) - 1


    # ----------------------------------------- Collocation Scheme ------------------------------------- #
    # read settings from config file
    N      = cfg["collocation"]["number_of_segments"]
    d      = cfg["collocation"]["num_points"]
    method = cfg["collocation"]["method"]
    tau_root, C, D, B = radau(d, method)

    # period of each mesh/segment
    mesh_T: float = (time_raw[-1] - time_raw[0]) / N

    # create a time vector that represents N + 1 points which represent the first
    # point of each mesh. Considering the continuity constraint, the first point of
    # segment/mesh n + 1 represents the last point of segment n.
    time_mesh = np.arange(time_raw[0], time_raw[-1] + dt_raw, mesh_T)

    # now, for each mesh, let's create a time vector that contains the times at the 3collocation points.
    time_col = np.zeros((N * d, ))
    time_col[::d] = time_mesh[:-1] + tau_root[1] * mesh_T      # 1st collocation point of every mesh
    time_col[1::d] = time_mesh[:-1] + tau_root[2] * mesh_T     # 2nd collocation point of every mesh
    time_col[2::d] = time_mesh[:-1] + tau_root[3] * mesh_T     # 3rd collocation point of every mesh



    # ------------ Ground Reaction Forces (if any) ----------- #
    if grf_path:
        grf_init = np.where((GRFs_raw['time'] < (time_raw[0] + dt_GRF / 2)) & (GRFs_raw['time'] >= (time_raw[0] - dt_GRF / 2)))[0][0]
        grf_end = np.where((GRFs_raw['time'] < (time_raw[-1] + dt_GRF / 2)) & (GRFs_raw['time'] >= (time_raw[-1] - dt_GRF / 2)))[0][0] + 1
        GRF_time = GRFs_raw['time'][grf_init:grf_end]
        GRFs_data = GRFs[grf_init:grf_end, grf_indices - 1]



    # ------------------------------------- Residuals & Torque Actuators --------------------------------------------- #
    # ----------- Residuals -------------- #
    res_sf = cfg["actuator_scale_factors"]["residuals"]                 # scale factor
    res_keys  = res_sf.keys()
    res_indices = []
    res_sf_value = {}
    res_sf_indices = {}
    for i, key in enumerate(list(res_keys)):
        items = res_sf[key]
        #
        res_sf_value[key] = items["magnitude"]
        res_sf_names = items["names"]
        res_sf_indices[key] = [i for i, name in enumerate(q_names) if name in res_sf_names]
        res_indices = res_indices + res_sf_indices[key]

    # ----------- Muscles -------------- #
    forceSet = model.getForceSet()
    muscleSet = forceSet.getMuscles()
    num_muscles = muscleSet.getSize()
    muscle_names: List[str] = get_item_names(muscleSet)
    MT_params = get_mt_parameters(model, muscle_names)


    muscle_dof_names: List[str] = cfg["dof_names_spanned_by_muscles"]
    num_muscle_dofs = len(muscle_dof_names)

    # find the indices of these DoFs wrt to the global index set of the generalised coordinates
    muscle_dof_indices = np.array([i for i, name in enumerate(q_names) if name in muscle_dof_names])


    # ----------- Torque Actuators -------------- #
    joint_moment_sf = cfg["actuator_scale_factors"]["joint_moments"]
    # we define one actuator per independent DoF here we check which coordinates the user has indicated as where to
    # assign the residual forces. We need to make sure that we correctly separate the residual forces from the actuator
    # forces/torques so that we can assign them to the correct set of constraints.
    exclude = set(res_indices) | set(no_tracking_indices) | set(muscle_dof_indices)
    actuator_indices = [i for i in q_ind_indices if i not in exclude]
    actuator_names = [name for i, name in enumerate(q_names) if i in actuator_indices]
    num_actuators = len(actuator_indices)

    # prepare actuator scale factor vector to be used later in the equality constraint set
    actuator_sf = np.ones((num_actuators, 1)) * joint_moment_sf

    # let's check whether the user has provided some actuators to be given a specific scale factor
    act_special_sf = cfg["actuator_scale_factors"]["special"]
    act_keys = act_special_sf.keys()
    for i, key in enumerate(list(act_keys)):
        items = act_special_sf[key]

        act_special_values = items["magnitude"]
        act_special_names = items["names"]

        act_special_indices = [i for i, name in enumerate(actuator_names) if name in act_special_names]
        actuator_sf[act_special_indices] = act_special_values


    # -------------------------------------- Bounds ------------------------------------------------------- #
    # ----------- X ----------- #
    special_q_bounds = cfg["bounds"]["special_qs"]
    x_bounds, x_scaling = set_x_bounds(Qs, time_raw, q_ind_indices, num_q_all, special_q_bounds)

    # ----------- torque actuators ----------- #
    torque_bounds, torque_scaling = set_torque_bounds(num_actuators)

    # ----------- residuals ----------- #
    residual_bounds = set_residual_bounds(cfg)

    # ----------- muscle state and controls ----------- #
    muscle_bounds, muscle_scaling = set_muscle_bounds(num_muscles)

    # -------------------------------------- Initial Guess ------------------------------------------------------- #
    x_guess = generate_x_guess(Qs, q_ind_indices, time_raw, time_mesh, time_col, num_q_all, x_scaling)
    torque_guess = generate_torque_guess(num_actuators, time_mesh, time_col)
    muscle_guess = generate_muscle_guess(num_muscles, time_mesh, time_col, muscle_scaling)


    # -------------------------------------- Experimental data scaling ----------------------------------------------- #
    # ---------- q ------------ #
    q_end_scaled_all = x_guess['x_end'][:, ::2]
    q_col_scaled_all = x_guess['x_col'][:, ::2]

    # retrieve just coordinates to track
    q_end_scaled = q_end_scaled_all[:, q_tracking_indices]
    q_col_scaled = q_col_scaled_all[:, q_tracking_indices]

    # ---------- q_dot ------------ #
    qdot_end_scaled_all = x_guess['x_end'][:, 1::2]
    qdot_col_scaled_all = x_guess['x_col'][:, 1::2]

    # retrieve just coordinates to track
    qdot_end_scaled = qdot_end_scaled_all[:, q_tracking_indices]
    qdot_col_scaled = qdot_col_scaled_all[:, q_tracking_indices]

    # -------------------------------------- NLP formulation ------------------------------------------------------- #
    opti = ca.Opti()

    #  Note: we're trying to optimise the state variable and the controls at the:
    #  1) collocation points
    #  2) mesh end-points: this is necessary to enforce continuity across segments
    #  We follow the same scheme for every subset of the design variable.

    # ----------- State ----------- #
    # interleaved [q, q_dot] for independent coords
    # 1) collocation points
    dims_x  = 2 * num_q_ind
    X_col = opti.variable(dims_x, N * d)
    opti.subject_to(x_bounds.X_lower <= (X_col <= x_bounds.X_upper))
    opti.set_initial(X_col, x_guess["x_col_ind"].T)

    # 2) mesh end-points
    X = opti.variable(dims_x, N + 1)
    opti.subject_to(x_bounds.X_lower <= (X <= x_bounds.X_upper))
    opti.set_initial(X, x_guess["x_end_ind"].T)

    # Muscle activation
    # 1) collocation points
    a_m_col = opti.variable(num_muscles, N * d)

    # 2) mesh end-points
    a_m = opti.variable(num_muscles, N + 1)


    # Torque actuators
    # 1) collocation points
    a_a_col = opti.variable(num_actuators, N * d)
    opti.subject_to(torque_bounds.a_a_lower <= (a_a_col <= torque_bounds.a_a_upper))
    opti.set_initial(a_a_col, torque_guess["a_a_col"].T)

    # 2) mesh end-points
    a_a = opti.variable(num_actuators, N + 1)
    opti.subject_to(torque_bounds.a_a_lower <= (a_a <= torque_bounds.a_a_upper))
    opti.set_initial(a_a, torque_guess["a_a"].T)

    print(f"States: {num_q_ind * 2 + num_actuators}")

    # -------------------- Controls ------------------------- #
    e_a = opti.variable(num_actuators, N + 1)
    opti.subject_to(torque_bounds.e_a_lower < (e_a <= torque_bounds.e_a_upper))
    opti.set_initial(e_a, torque_guess["e_a"].T)


    A_col = opti.variable(num_q_ind, N * d)
    opti.subject_to(x_bounds.Qsdotdot_lower_ind < (A_col <= x_bounds.Qsdotdot_upper_ind))
    opti.set_initial(A_col, x_guess["acc_col"].T)

    print(f"Controls: {num_actuators + num_q_ind}")


    func_map_in = [
        X[:, :-1],
        X_col,
        A_col,
        a_a[:, :-1],
        a_a_col,
        e_a[:, :-1],
        ca.MX(q_end_scaled[:-1, :].T),
        ca.MX(q_col_scaled.T),
        ca.MX(qdot_end_scaled[:-1, :].T),
        ca.MX(qdot_col_scaled.T)
    ]


    # scaled GRF data if the user provided it
    if grf_path:
        GRF_col = CubicSpline(GRF_time, GRFs_data)(time_col)

        scaling_GRF = np.maximum(np.abs(np.min(GRF_col, axis=0)), np.abs(np.max(GRF_col, axis=0))).reshape(1, -1)
        GRF_scaled = GRF_col / scaling_GRF

        func_map_in.append(ca.MX(GRF_scaled.T))


    # -------------------- Residuals ------------------------- #
    # we use the config file to register the residuals that the user specified
    residuals = cfg["bounds"]["residuals"]
    res_keys = residuals.keys()
    residual_vars = {}
    for key in list(res_keys):
        values = residuals[key]
        var_size: int = values[1]    # how many DoFs the residual forces are applied to

        # register decision variable
        residual_vars[key] = opti.variable(var_size, N * d)
        opti.subject_to(residual_bounds[key]['lower'] <= (residual_vars[key] <= residual_bounds[key]['upper']))
        opti.set_initial(residual_vars[key], np.zeros((var_size, N * d)))

        func_map_in.append(residual_vars[key])


    print(f"Decision variables: {opti.nx}")
    print(f"Bound constraints : {opti.ng}")


    # ------------------------------------------------------------------- #
    #The following section uses the "MX" object from the CasADi library
    #to register symbolic variables that will later be used in CasADi
    #"Function" objects, which dictate how these variables interact with
    #each other.

    # The suffix "k" indicates the variable at the mesh end-points.
    # "j", on the other hand, at the collocation points.
    # ------------------------------------------------------------------- #
    # ---------- state-vector ---------- #
    Xk = ca.MX.sym("Xk", 2 * num_q_ind)
    Xj = ca.MX.sym("Xj", 2 * num_q_ind, d)
    Xkj = ca.horzcat(Xk, Xj)

    # torque activation
    a_ak  = ca.MX.sym("a_ak",  num_actuators)
    a_aj  = ca.MX.sym("a_aj",  num_actuators, d)
    a_akj = ca.horzcat(a_ak, a_aj)


    # ---------- controls ---------- #
    # torque actuator excitation
    e_ak  = ca.MX.sym("e_ak", num_actuators)

    # generalised coordinate accelerations: remember, we're using implicit formulation. Therefore,
    #  accelerations are treated as "controls".
    Aj  = ca.MX.sym("Aj", num_q_ind, d)

    # ---------- experimental data to track ---------- #
    q_track_k  = ca.MX.sym("Qs_track_k",  num_q_tracking)
    q_track_j  = ca.MX.sym("Qs_track_j",  num_q_tracking, d)
    q_track_kj = ca.horzcat(q_track_k, q_track_j)

    qd_track_k  = ca.MX.sym("qd_track_k",  num_q_tracking)
    qd_track_j  = ca.MX.sym("qd_track_j",  num_q_tracking, d)
    qd_track_kj = ca.horzcat(qd_track_k, qd_track_j)

    func_in = [Xk, Xj, Aj, a_ak, a_aj, e_ak, q_track_k, q_track_j, qd_track_k, qd_track_j]

    if grf_path:
        num_grfs = GRF_col.shape[-1]
        GRF_track_j = ca.MX.sym("GRF_track_j", num_grfs, d)

        func_in.append(GRF_track_j)


    # ------------ Residuals ------------ #
    res_sym = {}
    for key in list(res_keys):

        values = residuals[key]
        var_size: int = values[1]    # how many DoFs the residual forces are applied to

        # register decision variable
        res_sym[key] = ca.MX.sym(f"{key}_j", var_size, d)
        func_in.append(res_sym[key])


    # Unscale X
    Xkj_nsc = ca.MX.zeros(*Xkj.shape)
    Xkj_nsc[0::2, :] = Xkj[0::2, :] * x_scaling.Qs[q_ind_indices]
    Xkj_nsc[1::2, :] = Xkj[1::2, :] * x_scaling.Qsdot[q_ind_indices]


    Aj_nsc  = Aj * x_scaling.Qsdotdot[q_ind_indices]


    # ---------------------- CasADi functions -------------------------- #
    # torque actuation dynamics
    actuator_dynamics_func = torque_activation_dynamics_casadi(num_actuators)

    # Cost-function: sum of squares
    J_q = sum_of_squares("q", num_q_tracking - num_special_indices)
    J_q_dot = sum_of_squares("q_dot", num_q_tracking)
    J_acc = sum_of_squares("acc", num_q_ind)

    #  check whether the user has defined any special q or q not to be tracked:
    #  we initialise this anyway so that even if the user hasn't provided
    #  any, we just add a 0 to the objective.
    q_no_tracking_term = ca.MX.zeros(1, 1)
    if no_tracking_names:
        J_no_tracking_term = sum_of_squares('no_tracking', num_no_tracking_indices)

    q_special_term = ca.MX.zeros(1, 1)
    if special_tracking_names:
        J_special = sum_of_squares('special', num_special_indices)

    # we do the same for the GRF term
    GRF_term = ca.MX.zeros(1, 1)
    if grf_path:
        J_GRF = sum_of_squares("GRF", num_grfs)


    # residual term
    J_res = {}
    res_term = {}
    for key in list(res_keys):
        values = residuals[key]
        var_size: int = values[1]

        # register residual cost functions
        J_res[key] = sum_of_squares(key, var_size)
        res_term[key] = ca.MX.zeros(1, 1)

    print(f"\nAll CasADi functions have been registered")

    #  Read the weights of the cost function from config file
    W = cfg["W"]

    #  ------------------------------------------------------------------- #
    # Prepare accumulators for what we want to numerically evaluate after optimisation:
    # cost function
    J  = ca.MX.zeros(1, 1)
    q_term = ca.MX.zeros(1, 1)
    q_dot_term = ca.MX.zeros(1, 1)
    acc_term = ca.MX.zeros(1, 1)
    residuals_term = ca.MX.zeros(1, 1)

    # initialise equality constraint list
    eq_constr  = []


    #  --------------------------------------------------------------------------------------------------------------- #
    # The following section is where the cost function and equality constraints (derivatives, rigid body dynamics,
    # torque actuator dynamics) are evaluated at the collocation points.
    # Here is where the CasADi builds the expression graph and the symbolic interdependency between the variables
    # that were previously registered.

    # The evaluation is done sequentially, over the 3 collocation points. It is written over the 3 collocation points
    # and later the CasADi parallelisation properties will be leveraged to apply the execution of this section of the
    # graph across multiple mesh segments at the same time, using multithreading.
    # ---------------------------------------------------------------------------------------------------------------- #
    for i in range(d):
        # current independent q and q_dot
        x_i  = Xkj_nsc[:, i + 1]

        # current independent accelerations
        acc_i = Aj_nsc[:, i]

        # apply kinematic couplings --> get full x and acc (independent and dependent)
        x_all_i, acc_all_i = apply_constraints(x_i, acc_i, q_names, q_ind_indices, path_to_config=kinematic_coupling_path)

        # evaluate external function: rigid body dynamics
        Ti = F(ca.vertcat(x_all_i, acc_all_i))

        # ------------------- Cost Function Terms ------------------- #
        # ---------- q ----------- #
        q_i = Xkj[::2, i + 1]

        if special_tracking_names:
            q_special_diff = q_i[special_q_indices_l] - q_track_kj[special_q_indices_l, i + 1]
            q_special_term =  q_special_term + W["q_special"] * B[i + 1] * J_special(q_special_diff) * mesh_T

        q_new_indices = [i for i in tracking_indices_l if i not in special_q_indices_l]
        q_diff = q_i[q_new_indices] - q_track_kj[indices_to_track, i + 1]
        q_term = q_term + W["q"] * B[i + 1] * J_q(q_diff) * mesh_T

        # regularisation term for q's not being tracked
        if no_tracking_names:
            q_no_tracking_term = q_no_tracking_term + W["q_no_tracking"] * B[i + 1] * J_no_tracking_term(q_i[no_tracking_indices_l]) * mesh_T


        # ---------- q dot ----------- #
        q_dot_i = Xkj[1::2, i + 1]
        q_dot_diff = q_dot_i[tracking_indices_l] - qd_track_kj[:, i + 1]
        q_dot_term = q_dot_term + W["q_dot"] * B[i + 1] * J_q_dot(q_dot_diff) * mesh_T

        # ---------- accelerations: regularisation term ----------- #
        acc_term = acc_term + W["acc"] * B[i + 1] * J_acc(Aj[:, i]) * mesh_T


        # ---------- GRFs (if any) ----------- #
        if grf_path:
            Ti_GRF_scaled = Ti[num_q_all + dll_grf_idx, 0] / scaling_GRF.T
            GRF_diff = Ti_GRF_scaled - GRF_track_j[:, i]
            GRF_term = GRF_term + W["GRF"] * B[i + 1] * J_GRF(GRF_diff) * mesh_T


        # ---------- residuals ----------- #
        for key in list(res_keys):
            J_func = J_res[key]
            res_term[key] = res_term[key] + W["residuals"] * B[i + 1] * J_func(res_sym[key][:, i]) * mesh_T

            residuals_term = residuals_term + res_term[key]


        # ---------- add them up ----------- #
        J = q_term + q_special_term + q_no_tracking_term + q_dot_term + acc_term + GRF_term + residuals_term


        # ------------------------------------------------------------------------------------------------------------ #
        #                                            Equality constraints                                              #
        # ------------------------------------------------------------------------------------------------------------ #
        # Rigid body dynamics:

        #  1) the derivative of the positional data at the current collocation point must be equal to the q_dot
        #  at the same collocation point and the derivative of 1_dot must be equal to the acceleration.
        #  To do this, we use the C matrix from the collocation scheme to compute the approximation derivatives.
        #  This is because, technically, IPOPT treats the two discretised trajectories (q, q_dot and acc), and therefore,
        #  the individual collocation points, independently.Therefore, we need to impose this equality
        #  constraint to make sure that the physics is correct.
        # approximate derivative of q
        qdot_nsc_approx = Xkj_nsc[::2, :] @ C[:, i + 1]       # [num_q, d + 1] @ [d+1, 1] -> [num_q, 1]
        # approximate derivative of q_dot
        acc_nsc_approx = Xkj_nsc[1::2, :] @ C[:, i + 1]

        # retrieve unscaled decision variables at the current collocation point
        qdot_nsc = Xkj_nsc[1::2, i + 1]

        eq_constr.append((mesh_T * qdot_nsc - qdot_nsc_approx) / x_scaling.Qs[q_ind_indices])
        eq_constr.append((mesh_T * acc_i - acc_nsc_approx) / x_scaling.Qsdot[q_ind_indices])


        # Torque activation dynamics:

        #  we impose the constraint that the derivative of the torque
        #  activation at the current collocation point must be equal to the
        #  analytical derivative that is computed via the activation
        #  dynamics formula.
        a_a_dot_approx = a_akj @ C[:, i + 1]          # [num_actuators, d + 1] @ [d+1, 1] -> [num_actuators, 1]

        # torque activation dynamics (explicit formulation)
        a_a_dot = actuator_dynamics_func(e_ak, a_akj[:, i + 1])
        eq_constr.append((mesh_T * a_a_dot - a_a_dot_approx))


        # Path constraints
        # --------------------------------------------------------------- %
        # here, we want to impose the constraint that
        # Computed torque from CasADi should be equal to the net moments
        # coming out of the OpenSim model. We do this only for the DoFs
        # that aren't spanned by the muscles.
        Ti_torques = Ti[actuator_indices, 0] / actuator_sf
        eq_constr.append(Ti_torques - a_akj[:, i + 1])

        # Also, we want to make sure that everything is dynamically consistent and therefore, we want to make sure that
        # the residuals computed from the current kineamatics is below the threshold that the user defined.
        for key in list(res_keys):
            res_tau = Ti[res_sf_indices[key], 0] / res_sf_value[key]
            eq_constr.append(res_sym[key][:, i] - res_tau)


    eq_constr  = ca.vertcat(*eq_constr)

    # Now we define a CasADi function that takes in the design variables at
    # the collocation points and outputs the cost function (J) and the sets
    # of constraints.
    func_out = [eq_constr, J, q_term, q_dot_term, acc_term]
    f_coll = ca.Function('f_coll', func_in, func_out)

    # save expression graph

    # register function as parallel form across the number of segments of
    # the trajectory.
    parallel_mode = "thread"
    num_threads   = 8
    f_coll_map = f_coll.map(N, parallel_mode, num_threads)

    # evaluate
    eq_constr_all, J_all, q_term_all, q_dot_term_all, acc_term_all = f_coll_map(*func_map_in)

    #  add constraints to opti struct
    opti.subject_to(eq_constr_all == 0)

    lbg = np.array(ca.evalf(opti.lbg)).flatten()
    ubg = np.array(ca.evalf(opti.ubg)).flatten()
    n_eq = np.sum(lbg == ubg)
    n_ineq = np.sum(lbg != ubg)

    print(f"Equality constraints before continuity: {n_eq}")
    print(f"Inequality constraints: {n_ineq}")
    print(f"Decision variables: {opti.nx}")


    # --------------------------------------------------------------------------------------------------------------- #
    #                                              equality constraints                                               #
    # --------------------------------------------------------------------------------------------------------------- #
    # because in direct collocation we treat each discrete point of the trajectory as a variable to optimise, we must
    # constrain how IPOPT can place these points in the y-axis. There is a significant risk of discontinuity in between
    # segments. Therefore, we must ensure that the segments along the trajectory are continous and the way to do it is
    # to impose as a constraint that the last point of a segment must be equal to the first point of its subsquent one.

    q_end = X[::2, :]
    q_col = X_col[::2, :]

    qdot_end = X[1::2, :]
    qdot_col = X_col[1::2, :]

    col_indices = np.array([0, 1, 2])
    # loop over the number of segments
    for k in range(N):
        q_kj = ca.horzcat(q_end[:, k], q_col[:, col_indices])
        qdot_kj = ca.horzcat(qdot_end[:, k], qdot_col[:, col_indices])
        act_kj = ca.horzcat(a_a[:, k], a_a_col[:, col_indices])

        opti.subject_to(q_end[:, k + 1] == q_kj @ D)
        opti.subject_to(qdot_end[:, k + 1] == qdot_kj @ D)
        opti.subject_to(a_a[:, k + 1] == act_kj @ D)

        col_indices += 3

    lbg = np.array(ca.evalf(opti.lbg)).flatten()
    ubg = np.array(ca.evalf(opti.ubg)).flatten()
    n_eq = np.sum(lbg == ubg)

    print(f"Final number of equality constraints: {n_eq}")

    # sum objective function over all mesh segments
    J_sum = ca.sum2(J_all)

    # --------------------------------------------------------------------------------------------------------------- #
    #                                                  NLP Solver                                                     #
    # --------------------------------------------------------------------------------------------------------------- #
    opti.minimize(J_sum)
    options = {}
    options['ipopt'] = {}
    options['ipopt']['hessian_approximation'] = 'limited-memory'
    options['ipopt']['mu_strategy'] = 'adaptive'
    options['ipopt']['max_iter'] = cfg["optimiser"]["max_iters"]
    tolerance = cfg["optimiser"]["tolerance"]
    options['ipopt']['tol'] = 1 * 10 ** (-tolerance)
    options['ipopt']['print_timing_statistics'] = 'yes'
    options['ipopt']['nlp_scaling_method'] = 'none'
    options['ipopt']['obj_scaling_factor'] = 1
    options['ipopt']['print_level'] = 5
    opti.solver('ipopt', options)

    # --------------------------------------------------------------- #
    # Solve problem
    w_opt, stats, g_opt, lambda_x, lambda_g = solve_NLP(opti, options)

    return w_opt, stats, g_opt, lambda_x, lambda_g

if __name__ == "__main__":
    model_path = r"C:\Users\ab3758\Documents\PhD\msk\P5\scaled_model_rot.osim"
    ik_path = r"C:\Users\ab3758\Documents\PhD\msk\P5\ik\P05R0002_ik_rot.mot"
    dll_path = r'C:\Users\ab3758\Documents\projects\tackling_msk\ocp\dlls\P05\build\RelWithDebInfo\P05.dll'
    config_filepath = r"C:\Users\ab3758\Documents\projects\tackling_msk\ocp\configs\tackling/config_rot.json"
    kinematic_coupling_path = r"C:\Users\ab3758\Documents\projects\tackling_msk\ocp\configs\tackling\kinematic_coupling_config.json"

    output_dir = r"C:\Users\ab3758\Documents\PhD\msk\P5\opt\py"
    grf_path = r"C:\Users\ab3758\Documents\PhD\msk\P5\grf\P05R0002.mot"

    muscle_driven(model_path, ik_path, dll_path, config_filepath, kinematic_coupling_path,output_dir, grf_path)