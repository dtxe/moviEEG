function [] = eegmov_readpointsfile(path)

fpoints = fopen(path);

data = textscan(fpoints, '%[ABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890 ]=%f,%f');
out.points.x = data{2};
out.points.y = data{3};
out.points.label = data{1};

out.num_channels = length(out.points.x);

%** Assign Channels **
out.num_samples = size(eegdata, 1);

% Match EEG channels from PAC to image coordinate labels
out.data = zeros(num_samples, out.num_channels);
for kk = 1:out.num_channels
    pos = strcmp(out.points.label(kk), ch_label);
    if max(pos) == 1
        out.data(:,kk) = eegdata(:, pos);
    else
        out.data(:,kk) = min(min(eegdata));
    end
end