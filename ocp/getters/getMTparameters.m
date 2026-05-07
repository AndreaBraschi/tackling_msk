% This function returns the muscle-tendon parameters of the muscles
% specified in muscleNames from the model in modelPath.
%
% Author: Antoine Falisse
% Date: 12/19/2018
% 
function params = getMTparameters(model, muscleNames)

import org.opensim.modeling.*

num_muscles = length(muscleNames);
params = zeros(5, num_muscles);

muscles = model.getMuscles();

for i = 1:num_muscles
   muscle = muscles.get(muscleNames{i});
   params(1,i) = muscle.getMaxIsometricForce();
   params(2,i) = muscle.getOptimalFiberLength();
   params(3,i) = muscle.getTendonSlackLength();
   params(4,i) = muscle.getPennationAngleAtOptimalFiberLength();
   params(5,i) = muscle.getMaxContractionVelocity()*params(2,i);
end

end