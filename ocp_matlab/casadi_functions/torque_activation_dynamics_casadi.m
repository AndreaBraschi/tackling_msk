function activation_dynamics_function = torque_activation_dynamics_casadi(muscle_model_path, num_actuators)
import casadi.*
addpath(muscle_model_path);

excitation  = SX.sym('excitation', num_actuators); % Excitation signals (input to the ODE)
activation  = SX.sym('activation', num_actuators); % Activations (ODE output)

dadt = torque_activation_dynamics(excitation, activation);

activation_dynamics_function = Function('torque_activation_dynamics_casadi', ...
    {excitation, activation}, {dadt});
end