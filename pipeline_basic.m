% Runs prospective cases through the eegmov pipeline
%
% Simeon Wong
% 2013 June 21

%% Import toolbox
addpath ./EEGmovie_pipeline
addpath ./brain-connectivity-toolbox

    
%% Parameters
FREQ_RANGE = [1 4; 4 7; 8 12; 13 32; 40 80; 80 150; 150 200;];
RESOLUTION_DIV = 3;     % divide the resolution for gridfitting. Higher number is faster. Lower number is higher quality.

WND_LENGTH = 4000;      % length of each window
WND_SHIFT  = 200;       % offset of subsequent window from previous

BASESAVE = '/data5/sample/save/directory/';

CONN_METRIC = 'plv';

USE_BASELINE = true;       % whether to baseline with interictal period


for s = 1:length(subj)
    % Reset everything, except for needed variables
    close all
    clearvars -except FREQ_RANGE RESOLUTION_DIV WND_LENGTH WND_SHIFT BASESAVE subj s CURRSCRIPT CONN_METRIC USE_BASELINE

    %% Computations
    fprintf('\n\n---\nRunning subject %s\n', subj{s});
    
    % Load subject data
    eegdata = eegmov_loadeegdata(SUBJECT_DATABASE, subj{s}, 'edf');
    
    % Filter with specified frequencies
    filtdata = eegmov_firfilthilbert(eegdata, FREQ_RANGE);
    
    % Calculate connectivity
    plvdata = eegmov_calcconnectivity(filtdata, WND_LENGTH, WND_SHIFT, CONN_METRIC);
    
    % Calculate subsequent graph properties
    graphmetrics = eegmov_graphmetrics(@eigenvector_centrality_und, plvdata);
    
    cmin = zeros(plvdata.num_freqs+1, 1);
    cmax = zeros(plvdata.num_freqs+1, 1);
    
    % If baselining with interictal
    if USE_BASELINE
        fprintf('Baselining...\n');

        intereegdata = eegmov_loadeegdata(SUBJECT_DATABASE, subj{s}, 'interedf');
        interfiltdata = eegmov_firfilthilbert(intereegdata, FREQ_RANGE);
        interplvdata = eegmov_calcconnectivity(interfiltdata, 0, 0, CONN_METRIC);
        intergraphmetrics = eegmov_graphmetrics(@eigenvector_centrality_und, interplvdata);
        
        for fq = 1:plvdata.num_freqs
            fqmean = mean(interplvdata.connectivity(:,:,fq), 2);
            for kk = 1:plvdata.num_windows
                plvdata.connectivity(:,kk,fq) = plvdata.connectivity(:,kk,fq) - fqmean;
            end
        end
        
        for fq = 1:graphmetrics.num_freqs
            for kk = 1:graphmetrics.num_windows
                graphmetrics.data(:,kk,fq) = graphmetrics.data(:,kk,fq) - ...
                    intergraphmetrics.data(:,1,fq);
            end
        end
        
        shownegative = true;

        cmaxval = max(abs(prctile(graphmetrics.data(:), 5)), abs(prctile(graphmetrics.data(:), 95)));
        cmin(1:end-1) = -1*cmaxval;
        cmax(1:end-1) = cmaxval;
    else
        shownegative = false;
        
        cmin(1:end-1) = prctile(graphmetrics.data(:), 40);
        cmax(1:end-1) = prctile(graphmetrics.data(:), 95);
    end

    % Calculate HFO amplitude data
    hfodata = eegmov_firfilthilbert(eegdata, [150 0], filtdata.num_samples_cut);     
    hfodata = eegmov_getHFO(hfodata, WND_LENGTH, WND_SHIFT);
    
    % Set colour threshold for hfodata overlay to 20th to 95th percentile
    cmin(end) = prctile(hfodata(:), 20);
    cmax(end) = prctile(hfodata(:), 95);
    
    % Take HFO amplitude data and add to existing graph theory summary statistic data
    % for each electrode
    graphmetrics.data = cat(3, graphmetrics.data, hfodata);
    graphmetrics.num_freqs = graphmetrics.num_freqs + 1;
    clear hfodata
    
    % Set connectivity in plvdata to none for HFO box, since HFO
    % is not a measure of connectivity
    plvdata.num_freqs = plvdata.num_freqs + 1;
    plvdata.connectivity = cat(3, plvdata.connectivity, zeros(plvdata.num_pairs, plvdata.num_windows) + NaN);

    % Interpolate electrode summary statistic data (HFO, graph properties)
    % to create surface overlay over entire image
    griddat.x = 1:RESOLUTION_DIV:size(eegdata.img, 2);
    griddat.y = 1:RESOLUTION_DIV:size(eegdata.img, 1);
    griddat.z = eegmov_gridfit(eegdata.points, graphmetrics, griddat);
    
    % Set things up for the movie generation script
    % This passes required information in the ampData structure
    ampData.WND_SHIFT = WND_SHIFT;
    ampData.WND_LENGTH = WND_LENGTH;
    ampData.num_samples_cut = filtdata.num_samples_cut;
    ampData.data = eegdata.ampElectrodes;
    ampData.labels = eegdata.ampElectrodeLabels;
    
    %% Make a movie
    mov_savepath = strcat(BASESAVE,'C',subj{s},'_w',num2str(WND_LENGTH),'-',num2str(WND_SHIFT),'_ch67.mp4');
    
    eegmov_rendermovie(mov_savepath, eegdata.img, eegdata.points,           ...
        round(eegdata.srate/WND_SHIFT), 'numpx', 4, 'numpy', 2,             ...
        'imoverlay', griddat, 'plvmatrix', plvdata,                         ...
        'infoText', {'','1-4Hz | 4-7Hz | 8-12Hz | 13-32Hz','40-80Hz | 80-150Hz | 150-200Hz | HFO Amplitude'},...
        'ampData', ampData, 'cmin', cmin, 'cmax', cmax, 'alpha', 0.2, 'minThreshold', 97,...
        'indicators', eegdata.indicators, 'showNegative', shownegative);
    
    % Notify user via pushover that subject has finished running
    alertstring = strcat('Movie saved: ', mov_savepath);
    fprintf(alertstring);

end