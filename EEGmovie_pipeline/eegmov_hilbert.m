function [ out ] = eegmov_hilbert( eegin, FREQ_RANGE, varargin )
% EEGMOV_HILBERT Filters data using FIR filter, then performs
% hilbert transform
%
% Simeon Wong
% 2013 June 5

num_freqs = size(FREQ_RANGE, 1);

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

