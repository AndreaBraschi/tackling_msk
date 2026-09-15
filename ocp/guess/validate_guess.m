function validate_guess(guess)
    fields = fieldnames(guess);
    for i = 1:numel(fields)
        field = fields{i};
        val = guess.(field);
        
        nan_idx = find(isnan(val(:)));
        inf_idx = find(isinf(val(:)));
        
        if ~isempty(nan_idx)
            error('[Initial Guess Validation] "%s" contains NaN at indices: %s', ...
                field, mat2str(nan_idx'));
        end
        
        if ~isempty(inf_idx)
            error('[Initial Guess Validation] "%s" contains Inf at indices: %s', ...
                field, mat2str(inf_idx'));
        end
    end
    fprintf('[Initial Guess Validation] All fields passed.\n');
end