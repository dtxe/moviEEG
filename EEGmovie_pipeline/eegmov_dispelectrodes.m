function [  ] = eegmov_dispelectrodes( eegdata )
%EEGMOV_DISPELECTRODES Summary of this function goes here
%   Detailed explanation goes here

f = figure;
image(eegdata.img);
hold on
scatter(eegdata.points.x, eegdata.points.y, 12^2, 'LineWidth', 2, 'MarkerEdgeColor', 'k');
axis equal;

end

