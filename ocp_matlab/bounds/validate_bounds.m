function validate_bounds(bounds)
    fields = fieldnames(bounds);
    for i = 1:numel(fields)
        field = fields{i};
        
        subfields = {'lower', 'upper'};
        for j = 1:numel(subfields)
            subfield = subfields{j};
            val = bounds.(field).(subfield);
            
            nan_idx = find(isnan(val(:)));
            inf_idx = find(isinf(val(:)));
            
            if ~isempty(nan_idx)
                error('[Bounds Validation] "%s.%s" contains NaN at indices: %s', ...
                    field, subfield, mat2str(nan_idx'));
            end
            
            if ~isempty(inf_idx)
                error('[Bounds Validation] "%s.%s" contains Inf at indices: %s', ...
                    field, subfield, mat2str(inf_idx'));
            end
        end
    end
    fprintf('[Bounds Validation] All fields passed.\n');
end