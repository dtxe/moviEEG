function [ out ] = eegmov_fir1filthilbert( eegin, FREQ_RANGE, varargin )
% EEGMOV_FIRFILTHILBERT Filters data using FIR filter, then performs
% hilbert transform
%
% Updated with the MATLAB recommended fir1 filter
%
% Simeon Wong
% 2015 Aug 11

num_freqs = size(FREQ_RANGE, 1);

time_filt = tic;

% Based upon code from eeglab

% ** Parameters **
trans = 0.3;   % fractional width of transition zones

nyq = eegin.srate / 2;

num_channels = eegin.num_channels;
filtdata = complex(zeros(eegin.num_samples, num_channels, num_freqs), 0);
srate = eegin.srate;

maxfiltorder = 0;

parfor fq = 1:num_freqs
    
    % ** Frequency Bins **
    locutoff = FREQ_RANGE(fq,1);
    hicutoff = FREQ_RANGE(fq,2);
    
    % ** Filter characteristics **
    
    if locutoff == 0        % lowpass filter
        transition = hicutoff * trans;
        freqs = [hicutoff, hicutoff + transition];
        fir_n = kaiserord(freqs, [1 0], [0.075 0.3], srate);
        
        filtwts2 = fir1(fir_n, hicutoff*2/srate, 'low');
    elseif hicutoff == 0    % highpass filter
        transition = locutoff * trans;
        freqs = [locutoff - transition, locutoff];
        fir_n = kaiserord(freqs, [0 1], [0.3 0.075], srate);
        
        filtwts2 = fir1(fir_n, locutoff*2/srate, 'bandpass');
    else
        % Estimate filter order
        transition = mean([locutoff, hicutoff]) * trans;
        freqs = [locutoff - transition, locutoff, hicutoff, hicutoff + transition];
        fir_n = kaiserord(freqs, [0 1 0], [0.3 0.075 0.3], srate);

        % Make filter coefficients
        filtwts2 = fir1(fir_n, [locutoff hicutoff]*2/srate, 'bandpass');
    end
    
    % ** Compute filtered and hilbert-ed signal **    
    filttemp = zeros(eegin.num_samples, eegin.num_channels);
    filttemp2 = zeros(eegin.num_samples, eegin.num_channels);
    for ch = 1:eegin.num_channels
        filttemp(:,ch) = filtfilt(filtwts2,1, eegin.data(:,ch) );
    end

    for ch = 1:eegin.num_channels
        filttemp2(:,ch) = hilbert(filttemp(:,ch));
    end
    
    maxfiltorder = max(maxfiltorder, fir_n);
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

