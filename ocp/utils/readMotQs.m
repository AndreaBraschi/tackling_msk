function Qs = readMotQs(filepath, cutoff)
% -------------------------------------------------------------------------
% readMotFile
%   This function reads a .mot file.

% INPUTs:
%   - mot_filepath (str): full path to the .mot file.
%   - filter (char): whether you want to filter the structure that gets
%   created from the .mot file data

% OUTPUT:
%   - Qs: filtered coordinates
% -------------------------------------------------------------------------
    arguments
        filepath
        cutoff = []   % default to empty if not passed
    end    
    

    Qs = importdata(filepath);

    Qs.time = Qs.data(:, 1);
    translational_tags = {'tx', 'ty', 'tz'};
    
    % When dealing with .mot files (IK results), we have the joint
    % coordinates in degrees. OpenSim/Simbody require radians for their
    % internal operations.
    % Let's look for the set of indices in the column space of the mot
    % structure that corresponds to translational DoFs
    mask = cellfun(@(h) any(contains(h, translational_tags)), Qs.colheaders);

    % now convert to radians
    Qs.data(:, ~mask) = Qs.data(:, ~mask) * (pi / 180);

    % Low-pass filter only if the user has passed a 'yes' as input
    if ~isempty(cutoff)
        order = 2;
        fs = 1 / mean(diff(Qs.data(:, 1)));
        [a, b] = butter(order/2, cutoff./(0.5*fs), 'low');
        Qs.allfilt = Qs.data;
        Qs.allfilt(:, 2:end) = filtfilt(a, b, Qs.allfilt(:, 2:end));
    
    else
        Qs.allfilt = Qs.data; % Store original data if no filtering is applied
    
    end


end
