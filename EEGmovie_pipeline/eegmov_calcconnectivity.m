function [ out ] = eegmov_calcconnectivity( filtered, WND_LENGTH, WND_SHIFT, metric )
%EEGMOV_CALCCONNECTIVITY Calculates connectivity between sensors.

if strcmp(metric, 'pli')
    calcconn = @(phasediff) abs(mean(sign(phasediff)));
elseif strcmp(metric, 'plv')
    calcconn = @(phasediff) abs(mean(exp(1i .* phasediff)));
else
    error('Invalid connectivity metric!');
end

time_plv = tic;

num_channels = filtered.num_channels;
num_freqs = filtered.num_freqs;

if WND_LENGTH ~= 0 && WND_SHIFT ~= 0
    num_windows = floor((filtered.num_samples - filtered.num_samples_cut*2 - WND_LENGTH) / WND_SHIFT)+1;
    fx_wnd_begin = NaN;
    fx_wnd_end = NaN;
else
    num_windows = 1;
    fx_wnd_begin = 1;
    fx_wnd_end = filtered.num_samples - filtered.num_samples_cut * 2;
end

num_pairs = num_channels * (num_channels - 1) / 2;
plvmatrix = zeros(num_pairs, num_windows, filtered.num_freqs);

for fq = 1:num_freqs
    
    parfiltdata = filtered.data(:,:,fq);
    parfor wnd = 1:num_windows
        if WND_LENGTH ~= 0 && WND_SHIFT ~= 0
            wnd_begin = (wnd - 1) * WND_SHIFT + 1;
            wnd_end = wnd_begin + WND_LENGTH - 1;
        else
            wnd_begin = fx_wnd_begin;
            wnd_end = fx_wnd_end;
        end
        
        parplvmatrix = zeros(num_pairs,1);

        pp = 0;
        for a = 1:num_channels-1
            for b = a+1:num_channels
                pp = pp + 1;
                phasediff = angle(parfiltdata(wnd_begin:wnd_end,a)) - angle(parfiltdata(wnd_begin:wnd_end,b));
                parplvmatrix(pp) = calcconn(phasediff);
            end
        end
        
        plvmatrix(:,wnd,fq) = parplvmatrix;
    end
    
end

time_plv = toc(time_plv);
fprintf('PLV time elapsed: %.1fs\n', time_plv);


out.connectivity = plvmatrix;
clear plvmatrix
out.num_pairs = num_pairs;
out.num_channels = num_channels;
out.num_windows = num_windows;
out.num_freqs = filtered.num_freqs;

end

