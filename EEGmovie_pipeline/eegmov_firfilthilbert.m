function [ out ] = eegmov_firfilthilbert( eegin, FREQ_RANGE, varargin )
% EEGMOV_FIRFILTHILBERT Filters data using FIR filter, then performs
% hilbert transform
%
% Simeon Wong
% 2013 June 5

num_freqs = size(FREQ_RANGE, 1);

time_filt = tic;

% Based upon code from eeglab

% ** Parameters **
minfac = 3;     % min frequency cycles in filter
trans = 0.10;   % fractional width of transition zones
mintap = 15;

nyq = eegin.srate / 2;

num_channels = eegin.num_channels;
filtdata = complex(zeros(eegin.num_samples, num_channels, num_freqs), 0);
srate = eegin.srate;

maxfiltorder = 0;

parfor fq = 1:num_freqs
    
    % ** Frequency Bins **
    locutoff = FREQ_RANGE(fq,1);
    hicutoff = FREQ_RANGE(fq,2);
    
    filtorder2 = minfac*fix(srate/locutoff);
    if filtorder2 < mintap
        filtorder2 = mintap;
    end
    
    % ** Filter characteristics **
    
    if locutoff == 0        % lowpass filter
        f2=[0 hicutoff/nyq (1+trans)*hicutoff/nyq 1];
        m =[1 1            0                      0];
    elseif hicutoff == 0    % highpass filter
        f2=[0 (1-trans)*locutoff/nyq locutoff/nyq 1];
        m =[0 0                      1            1]; 
    else
        f2=[0 (1-trans)*locutoff/nyq locutoff/nyq hicutoff/nyq (1+trans)*hicutoff/nyq 1];
        m =[0 0                      1            1            0                      0]; 
    end
    
    filtwts2 = firls(filtorder2,f2,m);
    
    % ** Compute filtered and hilbert-ed signal **    
    filttemp = zeros(eegin.num_samples, eegin.num_channels);
    filttemp2 = zeros(eegin.num_samples, eegin.num_channels);
    for ch = 1:eegin.num_channels
        filttemp(:,ch) = filtfilt(filtwts2,1, eegin.data(:,ch) );
    end

    for ch = 1:eegin.num_channels
        filttemp2(:,ch) = hilbert(filttemp(:,ch));
    end
    
    maxfiltorder = max(maxfiltorder, filtorder2);
    filtdata(:,:,fq) = filttemp2;
end

time_filt = toc(time_filt);
fprintf('Filter time elapsed: %.1fs\n', time_filt);

% ** Cut off filter junk data on both sides and write to output struct **
if ~isempty(varargin)
    num_samples_cut = varargin{1};
else
    num_samples_cut = ceil(maxfiltorder / 2);
end

out.data = filtdata(num_samples_cut:end-num_samples_cut, :, :);
clear filtdata

out.num_freqs = num_freqs;
out.num_samples = eegin.num_samples;
out.num_channels = eegin.num_channels;
out.num_samples_cut = num_samples_cut;
out.freq_range = FREQ_RANGE;

end

