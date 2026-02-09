function guess = getGuess(Qs, nq, N, nMuscles,joints,scaling, d, GRF, F1,Options,mode)
% --------------------------------------------------------------------------
% getGuess
%   This function provides an inital guess for the design variables that is
%   taken directly from the experimental data.

% INPUTs:
%   - Qs: coordinates coming from IK.
% 
%   - nq (integer): number of coordinates

%   - scaling -
% 
%   - d: degree of the interpolating polynomial of the collocation scheme
%
% OUTPUT:
%   - guess -
%   * initial guess values for all optimisation variables
% 
% Original author: Antoine Falisse
% Original date: 12/19/2018
%
% Last edit by: 
% Last edit date: 
% --------------------------------------------------------------------------