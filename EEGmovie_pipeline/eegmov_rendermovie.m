function [  ] = eegmov_rendermovie( savepath, img, points, fps, varargin)
%EEGMOV_RENDERMOVIE Summary of this function goes here
%   Detailed explanation goes here

%% Check inputs
p = inputParser;

addRequired(p, 'savepath');
addRequired(p, 'img');
addRequired(p, 'points');
addRequired(p, 'fps');
addParamValue(p, 'numpx', 4);
addParamValue(p, 'numpy', 2);
addParamValue(p, 'plotTitles', []);
addParamValue(p, 'infoText', {});

addParamValue(p, 'maxLineWidth', 4.0);
addParamValue(p, 'minLineWidth', 0.5);
addParamValue(p, 'minThreshold', 97);
addParamValue(p, 'maxThreshold', 100);

addParamValue(p, 'negMinThreshold', 0);   
addParamValue(p, 'negMaxThreshold', 3); 
addParamValue(p, 'showNegative', false);

addParamValue(p, 'plvmatrix', []);
addParamValue(p, 'imoverlay', []);

addParamValue(p, 'lineAlpha', 0.5);
addParamValue(p, 'alpha', 0.3);

addParamValue(p, 'cmin', []);
addParamValue(p, 'cmax', []);

addParamValue(p, 'spHeight', 255);
addParamValue(p, 'spWidth', 330);
addParamValue(p, 'spPadHeight', 10);
addParamValue(p, 'spPadWidth', 10);
addParamValue(p, 'ampHeight', 250);

addParamValue(p, 'ampData', []);

addParamValue(p, 'qpLevel', 3);

addParamValue(p, 'indicators', []);

parse(p, savepath, img, points, fps, varargin{:});
clear varargin

numpx = p.Results.numpx;
numpy = p.Results.numpy;

fps = p.Results.fps;

max_linewidth = p.Results.maxLineWidth;
min_linewidth = p.Results.minLineWidth;

save_path = p.Results.savepath;
[filea,fileb,~] = fileparts(save_path);

TRANSPARANCY = p.Results.alpha;
LINE_ALPHA = p.Results.lineAlpha;

cmin = p.Results.cmin;
cmax = p.Results.cmax;

spHeight = p.Results.spHeight;
spWidth = p.Results.spWidth;
spPadHeight = p.Results.spPadHeight;
spPadWidth = p.Results.spPadWidth;
ampHeight = p.Results.ampHeight;

ampData = p.Results.ampData;

infoText = {strrep(fileb, '_', '\_'), p.Results.infoText{:}};                                  %#ok<CCAT>

qp_level = p.Results.qpLevel;               % for video encoding

indicators = p.Results.indicators;


% Things to plot
plot_connectivity = ~isempty(p.Results.plvmatrix);
plot_overlay = ~isempty(p.Results.imoverlay);

if ~plot_connectivity && ~plot_overlay
    error('Nothing to plot!');
end

if plot_connectivity
    num_channels = p.Results.plvmatrix.num_channels;
    num_freqs = p.Results.plvmatrix.num_freqs;
    num_windows = p.Results.plvmatrix.num_windows;
else
    num_freqs = size(p.Results.imoverlay.z, 4);
    num_windows = size(p.Results.imoverlay.z, 3);
end

if plot_connectivity
    plvmatrix = p.Results.plvmatrix.connectivity;
    minThreshold = p.Results.minThreshold;
    maxThreshold = p.Results.maxThreshold;
    
    showNegative = p.Results.showNegative;
    
    if showNegative
        negMinThreshold = p.Results.negMinThreshold;
        negMaxThreshold = p.Results.negMaxThreshold;
    end
end
if plot_overlay
    grid_z = p.Results.imoverlay.z;
    grid_x = p.Results.imoverlay.x;
    grid_y = p.Results.imoverlay.y;
    
    if isempty(cmin) || isempty(cmax)
        cmin = 0;
        cmax = prctile(grid_z(:), 99);
    end
    
    if isscalar(cmin)
        cmin(1:num_freqs) = cmin;
        cmax(1:num_freqs) = cmax;
    end
end

clear p


%% Render movie

fprintf('Making pretty pictures...\n');
time_picture = tic;

imgsize = size(img);

spGetX = @(col) 20+col*(spWidth+spPadWidth);
spGetY = @(row) ampHeight+row*(spHeight+spPadHeight);

if plot_connectivity
    if isscalar(minThreshold)
        img_threshold = prctile(plvmatrix(:), minThreshold); 
        img_threshold = ones(num_windows, num_freqs) .* img_threshold;
    elseif isvector(minThreshold)     % assume across frequencies
        img_threshold = minThreshold;
        img_threshold = reshape(img_threshold, 1, []);
        img_threshold = repmat(img_threshold, num_windows, 1);
    elseif ismatrix(minThreshold)
        img_threshold = minThreshold;
    else
        error('minThreshold unsupported');
    end
    prc_max = prctile(plvmatrix(:), maxThreshold);
    
    if showNegative
        if isscalar(negMaxThreshold)
            neg_img_threshold = prctile(plvmatrix(:), negMaxThreshold); 
            neg_img_threshold = ones(num_windows, num_freqs) .* neg_img_threshold;
        elseif isvector(negMaxThreshold)     % assume across frequencies
            neg_img_threshold = negMaxThreshold;
            neg_img_threshold = reshape(neg_img_threshold, 1, []);
            neg_img_threshold = repmat(neg_img_threshold, num_windows, 1);
        elseif ismatrix(negMaxThreshold)
            neg_img_threshold = negMaxThreshold;
        else
            error('negMaxThreshold unsupported');
        end
    
        neg_prc_min = prctile(plvmatrix(:), negMinThreshold);
    end
end

if ~isempty(ampData)
    
    ampoffset = max(abs(ampData.data(:))) .* 0.9;

    ampSize = size(ampData.data);
    for nne = 1:ampSize(2)
        ampData.data(:,nne) = ampData.data(:,nne) + (ampoffset .* (nne - 1));
    end

end
    
temp_path = strcat(filea,'/',fileb,'.avi');
save_path = strcat(filea,'/',fileb,'.mp4');

[data_x, data_y] = meshgrid(1:imgsize(2), 1:imgsize(1));

vid = VideoWriter(temp_path);
vid.Quality = 100;
vid.FrameRate = fps;
open(vid);

f = figure;

% Determine figure size, based on subplot size and number
figSize = [(spWidth+spPadWidth)*numpx+25, (spHeight+spPadHeight)*numpy+ampHeight+10];

% Ensure figure size is even, otherwise, ffmpeg whines about it
figSize(1) = figSize(1) + (mod(figSize(1), 2) ~= 0);
figSize(2) = figSize(2) + (mod(figSize(2), 2) ~= 0);

textprogressbar('Progress: ');

set(f, 'Position', [10 10 figSize]);
set(f, 'Renderer', 'painters');
for kk = 1:num_windows
    figure(f);
    clf reset;
    
    for fq = 1:num_freqs        
        ax = subplot(numpy, numpx, fq);
        col = mod(fq-1, numpx);
        row = (numpy-1)-floor((fq-1) / numpx);
        
        set(ax, 'Units', 'pixels', 'Position', ...
            [spGetX(col), spGetY(row), spWidth, spHeight]);        
        
        image(img);
        hold on
        scatter(points.x, points.y, 12^2, 'LineWidth', 1.5, 'MarkerEdgeColor', 'k');
        
        if plot_overlay
            data_z = interp2(grid_x, grid_y, grid_z(:,:,kk,fq), data_x, data_y);

            himg = imagesc(data_z);
            set(himg, 'AlphaData', TRANSPARANCY);
            caxis([cmin(fq) cmax(fq)]);
        end

        if plot_connectivity
            pp = 0;
            for a = 1:num_channels-1
                for b = a+1:num_channels
                    pp = pp + 1;
                    if ~isnan(plvmatrix(pp,kk,fq)) && (plvmatrix(pp, kk, fq) >= img_threshold(kk,fq))
                        % draw line
                        linewidth = (plvmatrix(pp, kk, fq) - img_threshold(kk,fq)) / (prc_max - img_threshold(kk,fq)) * (max_linewidth - min_linewidth) + min_linewidth;
                        patchline([points.x(a),points.x(b)], [points.y(a),points.y(b)], 'LineWidth', linewidth, 'EdgeAlpha', LINE_ALPHA, 'EdgeColor', 'b');
                    end
                    
                    if showNegative && ~isnan(plvmatrix(pp,kk,fq)) && (plvmatrix(pp, kk, fq) <= neg_img_threshold(kk,fq))
                        % draw line
                        linewidth = (neg_img_threshold(kk,fq) - plvmatrix(pp, kk, fq)) / (neg_img_threshold(kk,fq) - neg_prc_min) * (max_linewidth - min_linewidth) + min_linewidth;
                        patchline([points.x(a),points.x(b)], [points.y(a),points.y(b)], 'LineWidth', linewidth, 'EdgeAlpha', LINE_ALPHA, 'EdgeColor', 'm');
                    end
                end
            end
        end
        
%         title(sprintf('%d-%dHz', FREQ_RANGE(1,fq), FREQ_RANGE(2,fq)));
        set(ax, 'XTick', []);
        set(ax, 'YTick', []);
        
        axis equal
    end
    
    colorbar('SouthOutside');
    set(ax, 'Units', 'pixels', 'Position', ...
            [spGetX(col), spGetY(row), spWidth, spHeight]);
    
    %% Render video info box
    infoPosition = [20+col*(spWidth+spPadWidth)+20, 25] ./ figSize;
    infoSize = [290, ampHeight-100] ./ figSize;
    annotation('textbox', [infoPosition, infoSize], 'String', infoText);
    
    %% Draw amplitude plot
    if(~isempty(ampData))
        ampPosition = [25, 25]./figSize;
        ampSize = [spGetX(col)-40, ampHeight - 50]./figSize;
        ax = subplot('Position', [ampPosition, ampSize]);
        
        plot(ampData.data);
        
        if any(strcmp('labels',fieldnames(ampData)))
            legend(ampData.labels{:});
        end
        
        yaxis = ylim();

        % Draw rectangle
        hP = patch([(kk-1)*ampData.WND_SHIFT+ampData.num_samples_cut, (kk-1)*ampData.WND_SHIFT+ampData.num_samples_cut, (kk-1)*ampData.WND_SHIFT+ampData.WND_LENGTH+ampData.num_samples_cut, (kk-1)*ampData.WND_SHIFT+ampData.WND_LENGTH+ampData.num_samples_cut],...
                   [yaxis(1), yaxis(2), yaxis(2), yaxis(1)], 'r');
        set(hP, 'FaceAlpha', 0.5);
        xlim([0 length(ampData.data)]);
        
        if ~isempty(indicators)
            for iind = 1:length(indicators)
                line([1 1] * indicators(iind), yaxis, 'LineStyle', '--', 'Color', 'k')
            end
        end
        
        set(ax, 'YTick', []);
    end
       
    refresh(f);
    writeVideo(vid, getframe(f));
%     writeVideo(vid, hardcopy(gcf, '-Dzbuffer', '-r0'));
    textprogressbar(kk/num_windows*100);
end
close(f);
close(vid);
textprogressbar('done');

time_picture = toc(time_picture);
fprintf('Picture rendering took %.1f seconds.\n\n', time_picture);

[status,output] = system(sprintf('ffmpeg -i %s -vcodec libx264 -qp %d -y %s', temp_path, qp_level, save_path));
if status == 0
    disp('Encoding success!');
    delete(temp_path);
else
    disp(output);
    disp('Error encoding video.');
end

end

