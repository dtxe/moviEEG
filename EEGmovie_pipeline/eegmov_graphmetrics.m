function [ out ] = eegmov_graphmetrics( graphfunction, plvdata )
%EEGMOV_GRAPHMETRICS Calculates graph theory metrics on data specified by
%                    function handle `graphfunction`

disp('Computing graph theory metrics...');
time_graph = tic;

plvmatrix = plvdata.connectivity;
num_windows = plvdata.num_windows;
num_freqs = plvdata.num_freqs;
num_channels = plvdata.num_channels;
clear plvdata

graphmetric = zeros(num_channels, num_windows, num_freqs);
parfor fq = 1:num_freqs
    for wnd = 1:num_windows
        connmat = squareform(plvmatrix(:,wnd,fq), 'tomatrix');
        graphmetric(:,wnd,fq) = graphfunction(connmat);                     %#ok<PFBNS>
    end
end

time_graph = toc(time_graph);
fprintf('Graph metrics took %.1f seconds.\n\n', time_graph);

out.num_channels = num_channels;
out.num_windows = num_windows;
out.num_freqs = num_freqs;
out.metric = func2str(graphfunction);
out.data = graphmetric;

end

