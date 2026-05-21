function GRF = readMotGrf(filepath, cutoff)
% -------------------------------------------------------------------------
% readMotGrf
%   This function reads a .mot file that contains the GRF and store them 
%   into a structure.

% INPUTs:
%   - filepath (str): full path to the .mot GRF file.
%   - cutoff (int): optional input that dictates the cutoff frequency of
%   the GRF data.

% OUTPUT:
%   - GRF: Ground Reaction Force structure
% -------------------------------------------------------------------------
    arguments
        filepath
        cutoff = []   % default to empty if not passed
    end    

    GRF = importdata(filepath);
    GRF.time = GRF.data(:, 1);

    % Low-pass filter only if the user has passed a 'yes' as input
    if ~isempty(cutoff)
        order = 2;
        fs = 1 / mean(diff(GRF.data(:, 1)));
        [a, b] = butter(order/2, cutoff./(0.5*fs), 'low');
        GRF_filt = filtfilt(a, b, GRF.data(:, 2:end));
        GRF.data(:, 2:end) = GRF_filt;
        
    end
end
