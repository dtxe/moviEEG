function [ out ] = eegmov_loadeegdata( case_num )
% LOADEEGDATA Loads eeg data using information from case parameters file
% and builds data struct for use in EEGmovie pipeline
%
% Simeon Wong
% 2013 June 5

    run(strcat(pwd, '/load_case/c', case_num, '.m'));

    fpoints = fopen(POINTS_PATH);
    
    out.edfpath = EDF_PATH;

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
    
    % Match EEG channels for amplitude plots
    if exist('ampElectrodes', 'var')
        out.ampElectrodes = zeros(num_samples, length(ampElectrodes));
        for kk = 1:length(ampElectrodes)
            pos = strcmp(ampElectrodes(kk), ch_label);
            if max(pos) == 1
                out.ampElectrodes(:,kk) = eegdata(:, pos);
            else
                error('Amplitude plot electrode not found');
            end
        end
    end

    out.img = imread(IMG_PATH);

    % scale point locations on image
    out.imgsize = size(out.img);
    out.points.x = (out.points.x * out.imgsize(2)/2) + (out.imgsize(2)/2);
    out.points.y = (out.points.y * out.imgsize(1)/2 * -1) + (out.imgsize(1)/2);
    out.srate = srate;
    
    if exist('INTERICTAL_FILE', 'var')
        out.interfn = INTERICTAL_FILE;
    else
        out.interfn = [];
    end
    out.ampElectrodeLabels = ampElectrodes;
    
    if exist('INDICATORS', 'var')
        out.indicators = INDICATORS;
    else
        out.indicators = srate * 60;
    end

end

