function force_equilibrium = buildForceEquilibriumFunc(muscle_model_path, num_muscles, MTparams)

import casadi.*
addpath(muscle_model_path);

% Function for Hill-equilibrium
FTtilde     = SX.sym('FTtilde', num_muscles); % Normalized tendon forces
a           = SX.sym('a', num_muscles); % Muscle activations
dFTtilde    = SX.sym('dFTtilde', num_muscles); % Time derivative tendon forces
lMT         = SX.sym('lMT', num_muscles); % Muscle-tendon lengths
vMT         = SX.sym('vMT', num_muscles); % Muscle-tendon velocities
Hilldiff    = SX(num_muscles, 1); % Hill-equilibrium
FT          = SX(num_muscles, 1); % Tendon forces
Fce         = SX(num_muscles, 1); % Contractile element forces
Fiso        = SX(num_muscles, 1); % Normalized forces from force-length curve
vMmax       = SX(num_muscles, 1); % Maximum contraction velocities

% Parameters of force-length-velocity curves
load Fvparam
load Fpparam
load Faparam

for m = 1:num_muscles
    [Hilldiff(m), FT(m), Fce(m), Fiso(m), vMmax(m)] = ...
        ForceEquilibrium_FtildeState_GC(a(m), FTtilde(m), dFTtilde(m),...
        lMT(m), vMT(m), MTparams(:,m), Fvparam, Fpparam, Faparam);
end

force_equilibrium = Function('force_equilibrium',{a, FTtilde, dFTtilde, lMT, vMT},...
    {Hilldiff, FT, Fce, Fiso, vMmax});
end