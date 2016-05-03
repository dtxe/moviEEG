function [ eegdata, EMD_allfiles ] = eegmov_getEMD( subj, varargin )
% EEGMOV_GETEMD Get the empirical mode decomposition for the given subject
%   If the EMD data already exists and is precomputed, use that data
%   instead.

%% Parse parameters
p = inputParser;
addParamValue(p, 'dir', []);
addParamValue(p, 'verbose', false);
addParamValue(p, 'overwrite', false);

parse(p, varargin{:});

%% Calculations

eegdata = eegmov_loadeegdata(subj);

if ~isempty(p.Results.dir)
    if p.Results.dir(end) ~= '/'
        basepath = strcat(p.Results.dir, '/');
    else
        basepath = p.Results.dir;
    end
    
    emddatapath = strcat(basepath, 'emd_', subj, '.mat');
    
    if exist(emddatapath, 'file') && ~p.Results.overwrite
        load(emddatapath, 'EMD_allfiles');
        
        if exist('EMD_allfiles', 'var')
            if p.Results.verbose
                fprintf('Precomputed EMD files found. Using precomputed values instead.\n');
            end
            
            return
        end
    end
end

% If EMD doesn't exist, recalculate!
eegsig = eegdata.data;

EMD_allfiles = {};
parfor kk = 1:eegdata.num_channels
    if p.Results.verbose
        fprintf('Running channel %d\n', kk);
    end

    EMD_allfiles{kk} = rParabEmd_nooutput(eegsig(:,kk),  40,40,1);        % Calculate EMD
end

if ~isempty(p.Results.dir)
    save(emddatapath, 'EMD_allfiles', '-v7.3');
end

end