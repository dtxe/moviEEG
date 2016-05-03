function [ hfodata ] = eegmov_getHFO( filtered, WND_LENGTH, WND_SHIFT )
%EEGMOV_GETHFO Windows HFO filtered data, and takes mean amplitude

num_windows = floor((filtered.num_samples - filtered.num_samples_cut*2 - WND_LENGTH) / WND_SHIFT)+1;

hfodata = zeros(filtered.num_channels, num_windows);
for wnd = 1:num_windows
    wnd_begin = (wnd - 1) * WND_SHIFT + 1;
    wnd_end = wnd_begin + WND_LENGTH - 1;
    
    hfodata(:,wnd) = mean(abs(filtered.data(wnd_begin:wnd_end,:,1)), 1);
end

end

