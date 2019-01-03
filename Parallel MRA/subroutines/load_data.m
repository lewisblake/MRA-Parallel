function [ data, regressionModel, domainBoundaries, predictionVector, theta, varEps ] = load_data(dataSource, nXGrid, nYGrid, offsetPercentage, verbose)
%% LOAD_DATA Loads data from various sources
%   Data files are loaded as a function of a variable called dataType
%
%   Input: dataSource (char), nXGrid (double), nYGrid (double),
%   offsetPercentage (double), verbose (boolean)
%
%   Output: data (matrix), regressionModel (LinearModel), domainBoundaries
%   (vector), predictionVector (matrix), theta (vector), varEps (double)
%%
switch dataSource   
    case 'satellite'
        %% User Input
        % Change as needed
        load('./Data/satelliteData.mat')
        
        % Values of parameters of covariance function
        theta = [5.57,0.12]; varEps = 0.01;
        
    case 'simulated'
        %% User Input
        % Change as needed
        load('./Data/simulatedData.mat')
        
        % Values of parameters of covariance function
        theta = [8.13,0.72]; varEps = 0.1;
        
    otherwise 
        error('Error. Specified dataType is not a valid data set.');
end
if verbose
    disp('Loading data complete');
end
% Determine the boundaries of the domain spanded by the data.
xMin0 = min(lon);
xMax0 = max(lon);
yMin0 = min(lat);
yMax0 = max(lat);
domainBoundaries = [xMin0, xMax0, yMin0, yMax0];

% Make prediction grid
if nXGrid && nYGrid > 0 % If user defines a prediction grid
    xPredictionVec = linspace(xMin0 + offsetPercentage*(xMax0 - xMin0), xMax0 - offsetPercentage*(xMax0 - xMin0), nXGrid);
    yPredictionVec = linspace(yMin0 + offsetPercentage*(yMax0 - yMin0), yMax0 - offsetPercentage*(yMax0 - yMin0), nYGrid);
    [xPredGridLocs, yPredGridLocs] = meshgrid(xPredictionVec, yPredictionVec);
else
    warning('Prediction grid set to be empty.')
    xPredGridLocs = [];
    yPredGridLocs = [];
end

% Find observation locations.
logicalInd = ~isnan(obs);
% Declare predictionVector
predictionVector = [xPredGridLocs(:),yPredGridLocs(:)];

% Assign lon, lat and observations to data matrix.
data(:,1) = lon(logicalInd);
data(:,2) = lat(logicalInd);
data(:,4) = obs(logicalInd);

% Detrend data.
regressionModel = fitlm(data(:,1:2),data(:,4), 'linear');
residuals = table2array(regressionModel.Residuals(:, 1));
data(:,3) = residuals;
end