function bounds = set_residualBounds(config_struct)

% This function assign the user defined residual forces.

%
% Inputs:
%   - config_struct (struct): 


residuals = config_struct.bounds.residuals;
keys      = fieldnames(residuals);

for i = 1:length(keys)
    key   = keys{i};
    values = residuals.(key);
    residual_value = values(1);   % magnitude of residual force
    num_coords = values(2);       % how many DoFs it is applied to
    
    lower = -residual_value * ones(1, num_coords);
    upper = residual_value * ones(1, num_coords);
  
    
    scaling  = max(abs(lower), abs(upper)); 
    bounds.(key).lower = (lower)./scaling;
    bounds.(key).upper = (upper)./scaling;
        

end


% -------------------------- Validate struct -------------------------- %
validate_bounds(bounds);


end