function [] = create_plots(data, predictionLocations, predictionMean, predictionVariance, verbose, displayPlots, savePlots, plotsFilePath)
%% CREATE_PLOTS creates the plots for the 'prediction' calculationType when plotting is true
%
%   Input: data, predictionLocations, predictionMean, predictionVariance,
% verbose, displayPlots, savePlots, plotsFilePath
%
%   Output: [] (empty vector)
%

if verbose
    disp('Beginning plotting')
end

% Figure 1
figure
if ~displayPlots
    set(gcf,'Visible', 'off');
end
scatter(data(:,1), data(:,2), 5, data(:,4),'square', 'filled');
colormap(parula)
colorbar
[cmin,cmax] = caxis;
caxis([cmin, cmax])
title('Observations')
xlabel('x coordinates')
ylabel('y coordinates')
if savePlots
    saveas(gcf, fullfile(plotsFilePath, 'observations'), 'png');
end

% Figure 2
figure
if ~displayPlots
    set(gcf,'Visible', 'off');
end
scatter(predictionLocations(:,1), predictionLocations(:,2), 5, predictionMean, 'square', 'filled');
colormap(parula)
colorbar
caxis([cmin, cmax])
title('Prediction mean')
xlabel('x coordinates')
ylabel('y coordinates')
if savePlots
    saveas(gcf, fullfile(plotsFilePath, 'prediction_mean'), 'png');
end

% Figure 3
figure
if ~displayPlots
    set(gcf,'Visible', 'off');
end
scatter(predictionLocations(:,1), predictionLocations(:,2), 5, predictionVariance, 'square', 'filled');
colormap(flip(autumn))
colorbar
title('Prediction variance')
xlabel('x coordinates')
ylabel('y coordinates')
if savePlots
    saveas(gcf, fullfile(plotsFilePath, 'prediction_variance'), 'png');
end

if verbose % Display progress indicators
    disp('Plotting completed')
end
end

