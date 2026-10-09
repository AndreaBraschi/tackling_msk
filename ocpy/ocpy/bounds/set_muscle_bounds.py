import numpy as np

def set_muscle_bounds(num_muscles: int):

    scaling = {}
    bounds = {}

    # ----------- activation ----------- #
    a_lower = np.zeros((num_muscles, 1))       # activation can't be negative
    a_upper = np.ones((num_muscles, 1))
    bounds['a_lower'] = a_lower
    bounds['a_upper'] = a_upper

    # ----------- activation time derivative ----------- #
    t_act = 0.015
    t_deact = 0.06
    va_lower = - ((1 / 100) * np.ones((num_muscles, 1))) / (np.ones((num_muscles, 1)) * t_deact)
    va_upper = ((1 / 100) * np.ones((num_muscles, 1))) / (np.ones((num_muscles, 1)) * t_act)
    bounds['va_lower'] = va_lower
    bounds['va_upper'] = va_upper

    scaling['va'] = 100

    # ----------------- Muscle Tendon Unit (MTU) Force ------------ #
    mtu_f_lower = np.zeros((num_muscles, 1))        # scalar force, can't be negative. Negative force is determined by direction
    mtu_f_upper = 5 * np.ones((num_muscles, 1))

    scaling['mtu_f'] = np.maximum(np.abs(mtu_f_lower), np.abs(mtu_f_upper))
    bounds['mtu_f_lower'] = mtu_f_lower / scaling['mtu_f']
    bounds['mtu_f_upper'] = mtu_f_upper / scaling['mtu_f']


    # ----------------- MTU force time derivative --------------- #
    mtu_df_lower = -np.ones((num_muscles, 1))
    mtu_df_upper = np.ones((num_muscles, 1))
    bounds['mtu_df_lower'] = mtu_df_lower
    bounds['mtu_df_upper'] = mtu_df_upper

    scaling['mtu_df'] = 100

    return bounds, scaling

