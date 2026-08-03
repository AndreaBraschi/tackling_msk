function [GRFs_num, pelvis_res_num, eq_constr_num, mt_ft_num, mt_dm_num] =  evaluate_guess(opti, J_all, q_term_all, q_dot_term_all, q_dot_dot_term_all, ...
   GRF_term_all, GRM_term_all, act_term_all, act_der_term_all, ...
   dMTf_dt_term_all, pelvis_term_all, GRF_eval_all, pelvis_res_eval_all, eq_constr_all, mt_ft_all, mt_dm_all)

    fprintf('q_term        : %f\n', opti.debug.value(sum(q_term_all), opti.initial));
    fprintf('q_dot_term    : %f\n', opti.debug.value(sum(q_dot_term_all), opti.initial));
    fprintf('q_dot_dot_term: %f\n', opti.debug.value(sum(q_dot_dot_term_all), opti.initial));
    fprintf('GRF_term      : %f\n', opti.debug.value(sum(GRF_term_all), opti.initial));
    fprintf('GRM_term      : %f\n', opti.debug.value(sum(GRM_term_all), opti.initial));
    fprintf('act_term      : %f\n', opti.debug.value(sum(act_term_all), opti.initial));
    fprintf('act_der_term  : %f\n', opti.debug.value(sum(act_der_term_all), opti.initial));
    fprintf('dMTf_dt_term  : %f\n', opti.debug.value(sum(dMTf_dt_term_all), opti.initial));
    fprintf('pelvis_term   : %f\n', opti.debug.value(sum(pelvis_term_all), opti.initial));
    fprintf('─────────────────────────────\n');
    fprintf('J total       : %f\n', opti.debug.value(sum(J_all), opti.initial));

    GRFs_num = opti.debug.value(GRF_eval_all, opti.initial);
    pelvis_res_num = opti.debug.value(pelvis_res_eval_all, opti.initial);
    eq_constr_num = opti.debug.value(eq_constr_all, opti.initial);
    mt_ft_num = opti.debug.value(mt_ft_all, opti.initial);
    mt_dm_num = opti.debug.value(mt_dm_all, opti.initial);


    fprintf('\n----------Equality constraints-----------\n');
    fprintf('velocity term        : %f\n', sum(eq_constr_num(1, :)));
    fprintf('acceleration term    : %f\n', sum(eq_constr_num(2, :)));
    fprintf('joint moment term: %f\n', sum(eq_constr_num(3, :)));
    fprintf('torque actuators activation term      : %f\n', sum(eq_constr_num(4, :)));
    fprintf('pelvis residuals term      : %f\n', sum(eq_constr_num(5, :)));
    fprintf('derivative of muscle activation term      : %f\n', sum(eq_constr_num(6, :)));
    fprintf('derivative of MT unit force term      : %f\n', sum(eq_constr_num(7, :)));
    fprintf('Hill error term      : %f\n', sum(eq_constr_num(8, :)));
    fprintf('─────────────────────────────\n');

end