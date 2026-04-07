function Qs = readMotQs(mot_filepath)
% -------------------------------------------------------------------------
% readMotFile
%   This function reads a .mot file.

% INPUTs:
%   - mot_filepath (str): full path to the .mot file.

% OUTPUT:
%   - Qs: filtered coordinates
% -------------------------------------------------------------------------
    mot = importdata(mot_filepath);

    Qs.time = mot.data(:,strcmp(mot.colheaders, {'time'}));
    Qs.all(:, 1) = Qs.time;
    Qs.colheaders{1,1} = 'time';

    count = 1;
    translational_tags = {"tx", "ty", "tz"};
    % loop through each column of the data
    for i = 1:size(mot.data, 2)
    count = count + 1;
    if length(split(mot.colheaders{i})) > 1 || any(ismember(split(mot.colheaders{i}), translational_tags))
        Qs.(mot.colheaders{i}) = Qsall.data(:,strcmp(Qsall.colheaders, mot.colheaders{i}));
    else
        Qs.(mot.colheaders{i}) = Qsall.data(:,strcmp(Qsall.colheaders, mot.colheaders{i})).*(pi/180);
    end
    Qs.all(:,count) = Qs.(mot.colheaders{i});
    Qs.colheaders(1,count) = Qsall.colheaders(1,strcmp(Qsall.colheaders, mot.colheaders{i}));

    end

    % Low-pass filter
    order = 2;
    cutoff_low = 15;
    fs=1/mean(diff(Qs.all(:,1)));
    [af,bf] = butter(order/2,cutoff_low./(0.5*fs),'low');
    Qs.allfilt = Qs.all;
    Qs.allfilt(:,2:end) = filtfilt(af,bf,Qs.allfilt(:,2:end));

end
