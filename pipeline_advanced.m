% Runs prospective cases through the eegmov pipeline
%
% Simeon Wong
% 2013 June 21

%% Import toolbox
if exist('eegmov_loadeegdata', 'file') ~= 2
    addpath /sdata3/simeon/prod/eegmovie/EEGmovie_pipeline
end

% Import runnotify
addpath /sdata3/simeon/prod/runnotify

%% Reset
clearvars -except subj
if matlabpool('size') ~= 0
    matlabpool('close');
end

% Get current script name for pipeline variants
[~,CURRSCRIPT,~] = fileparts(mfilename('fullpath'));

%% Parameters
FREQ_RANGE = [1 4; 4 7; 8 12; 13 32; 40 80; 80 150; 150 200;];
RESOLUTION_DIV = 3;

WND_LENGTH = 4000;      % length of each window
WND_SHIFT  = 200;       % offset of subsequent window from previous

BASESAVE = '/sdata2/simeon/test_data/prospective/';

if ~isempty(strfind(CURRSCRIPT, 'pli'))
    CONN_METRIC = 'pli';
else
    CONN_METRIC = 'plv';
end

% subj = {'120a', '120b', '120c'};


for s = 1:length(subj)
    close all
    clearvars -except FREQ_RANGE RESOLUTION_DIV WND_LENGTH WND_SHIFT BASESAVE subj s CURRSCRIPT CONN_METRIC
    matlabpool_custom

    %% Computations
    fprintf('\n\n---\nRunning subject %s\n', subj{s});
    eegdata = eegmov_loadeegdata(subj{s});
    filtdata = eegmov_firfilthilbert(eegdata, FREQ_RANGE);    
    plvdata = eegmov_calcconnectivity(filtdata, WND_LENGTH, WND_SHIFT, CONN_METRIC);
    graphmetrics = eegmov_graphmetrics(@eigenvector_centrality_und, plvdata);
    
    cmin = zeros(plvdata.num_freqs+1, 1);
    cmax = zeros(plvdata.num_freqs+1, 1);
    
    % If baselining with interictal
    if ~isempty(strfind(CURRSCRIPT, 'baseline'))
        fprintf('Baselining...\n');
        
        % Get interictal data
        if ~isempty(eegdata.interfn)
            interfn = eegdata.interfn;
        else
            interfn = strcat(subj{s}, 'inter');
        end
        
        intereegdata = eegmov_loadeegdata(interfn);
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

    
    hfodata = eegmov_firfilthilbert(eegdata, [150 0], filtdata.num_samples_cut);     
    hfodata = eegmov_getHFO(hfodata, WND_LENGTH, WND_SHIFT);
    
    cmin(end) = prctile(hfodata(:), 20);
    cmax(end) = prctile(hfodata(:), 95);
    
    graphmetrics.data = cat(3, graphmetrics.data, hfodata);
    clear hfodata
    
    graphmetrics.num_freqs = graphmetrics.num_freqs + 1;
    plvdata.num_freqs = plvdata.num_freqs + 1;
    plvdata.connectivity = cat(3, plvdata.connectivity, zeros(plvdata.num_pairs, plvdata.num_windows) + NaN);
    
    if ~isempty(strfind(CURRSCRIPT, 'dthresh'))
        % Get dynamic thresholds for connectivity plots
        imgthreshold = zeros(plvdata.num_windows, plvdata.num_freqs);
        for ff = 1:7
            for kk = 1:plvdata.num_windows
                imgthreshold(kk,ff) = prctile(plvdata.connectivity(:,kk,ff), 98.5);
            end
        end
    else
        imgthreshold = 95;
    end
    
    ampData.WND_SHIFT = WND_SHIFT;
    ampData.WND_LENGTH = WND_LENGTH;
    ampData.num_samples_cut = filtdata.num_samples_cut;
    ampData.data = eegdata.ampElectrodes;
    ampData.labels = eegdata.ampElectrodeLabels;

    griddat.x = 1:RESOLUTION_DIV:size(eegdata.img, 2);
    griddat.y = 1:RESOLUTION_DIV:size(eegdata.img, 1);
    griddat.z = eegmov_gridfit(eegdata.points, graphmetrics, griddat);

    matlabpool close
    
    %% Make a movie
    mov_savepath = strcat(BASESAVE,'C',subj{s},'_w',num2str(WND_LENGTH),'-',num2str(WND_SHIFT),'_',CURRSCRIPT,'_ch67.mp4');
    
    save(strcat('/data5/simeon/precomputed_eegmov/eegmov_', CURRSCRIPT, '_', subj{s}, '.mat'), ...
        'mov_savepath', 'WND_LENGTH', 'WND_SHIFT', 'eegdata', 'griddat', 'plvdata', 'ampData', 'cmin', 'cmax', ...
        'imgthreshold', 'graphmetrics', 'shownegative', '-v7.3');
    
    eegmov_rendermovie(mov_savepath, eegdata.img, eegdata.points,           ...
        round(eegdata.srate/WND_SHIFT), 'numpx', 4, 'numpy', 2,             ...
        'imoverlay', griddat, 'plvmatrix', plvdata,                         ...
        'infoText', {'','1-4Hz | 4-7Hz | 8-12Hz | 13-32Hz','40-80Hz | 80-150Hz | 150-200Hz | HFO Amplitude'},...
        'ampData', ampData, 'cmin', cmin, 'cmax', cmax, 'alpha', 0.2, 'minThreshold', imgthreshold,...
        'indicators', eegdata.indicators, 'showNegative', shownegative, 'isolateCh', 46);
    
    alertstring = strcat('Movie saved: ', mov_savepath);
    notify(alertstring);

end