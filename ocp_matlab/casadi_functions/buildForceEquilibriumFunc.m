function force_equilibrium_sym = buildForceEquilibriumFunc(pathMuscleModel, num_muscles, MTparameters)
% Function for Hill-equilibrium
import casadi.*
addpath(pathMuscleModel);

% read config
json_str = fileread(fullfile(pathMuscleModel, "config.json"));
config_struct = jsondecode(json_str);

% retrieve parameters of force-length-velocity curves
Fvparam = config_struct.("force_velocity_length").("Fvparam");
Fpparam = config_struct.("force_velocity_length").("Fpparam");
Faparam = config_struct.("force_velocity_length").("Faparam");

% inputs: those variables that will be numerically filled in.
FTtilde     = SX.sym('FTtilde', num_muscles); % Normalized tendon forces
a           = SX.sym('a', num_muscles); % Muscle activations
dFTtilde    = SX.sym('dFTtilde', num_muscles); % Time derivative tendon forces
lMT         = SX.sym('lMT', num_muscles); % Muscle-tendon lengths
vMT         = SX.sym('vMT', num_muscles); % Muscle-tendon velocities

% outputs: pre-register dimensionality. These will be the variables that
% will determine how the inputs are related to one another.
Hilldiff    = SX(num_muscles, 1); % Hill-equilibrium
FT          = SX(num_muscles, 1); % Tendon forces
Fce         = SX(num_muscles, 1); % Contractile element forces
Fiso        = SX(num_muscles, 1); % Normalized forces from force-length curve
vMmax       = SX(num_muscles, 1); % Maximum contraction velocities

for m = 1:num_muscles
    [Hilldiff(m), FT(m), Fce(m), Fiso(m), vMmax(m)] = ...
        force_equilibrium(a(m), FTtilde(m), dFTtilde(m),...
        lMT(m), vMT(m), MTparameters(:,m), Fvparam, Fpparam, Faparam);
end

force_equilibrium_sym = Function('force_equilibrium', {a, FTtilde, dFTtilde, lMT, vMT},...
    {Hilldiff, FT, Fce, Fiso, vMmax});
end