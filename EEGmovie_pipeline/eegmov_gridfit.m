function [ grid_z ] = eegmov_gridfit( points, surfdata, grid, varargin )
%EEGMOV_GRIDFIT Fit surface to given z-values of electrode grid

%% Parse inputs
p = inputParser;

addRequired(p, 'points');
addRequired(p, 'surfdata');
addRequired(p, 'grid');
addOptional(p, 'REGULARIZER_MODE', 'gradient');

parse(p, points, surfdata, grid, varargin);

points = p.Results.points;
surfdata = p.Results.surfdata;
grid = p.Results.grid;
REGULARIZER_MODE = p.Results.REGULARIZER_MODE;
clear p;

%% Calculate
grid_z = zeros(length(grid.y), length(grid.x), surfdata.num_windows, surfdata.num_freqs);

disp('Fitting grid...');
time_gridfit = tic;

for f=1:surfdata.num_freqs
    freqdata = surfdata.data(:,:,f);
    parfor kk=1:surfdata.num_windows
        [grid_z(:,:,kk,f), ~, ~] = gridfit(points.x, points.y, freqdata(:,kk), grid.x, grid.y, 'regularizer', REGULARIZER_MODE);
    end
end

time_gridfit = toc(time_gridfit);
fprintf('Gridfitting took %.1f seconds.\n\n', time_gridfit);

end

