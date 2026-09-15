function validate_array(val, name)
    nan_idx = find(isnan(val(:)));
    inf_idx = find(isinf(val(:)));
    
    if ~isempty(nan_idx)
        error('[Array Validation] "%s" contains NaN at indices: %s', ...
            name, mat2str(nan_idx'));
    end
    
    if ~isempty(inf_idx)
        error('[Array Validation] "%s" contains Inf at indices: %s', ...
            name, mat2str(inf_idx'));
    end
end
