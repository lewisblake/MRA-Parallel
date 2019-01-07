%% Multi-resolution Approximation (MRA) main script
% Executes the MRA model for three calculationType options
%
% Authors: Lewis Blake, Colorado School of Mines (lblake@mines.edu)
%          Dorit Hammerling, Colorado School of Mines (hammerling@mines.edu) 
%
%% Load User Input
% Run user_input.m script to get variables into workspace
user_input;
%% Check User Input is valid
validate_user_input(calculationType, NUM_LEVELS_M, NUM_PARTITIONS_J, NUM_KNOTS_r, offsetPercentage, NUM_WORKERS, nXGrid, nYGrid, displayPlots, savePlots, verbose, resultsFilePath, plotsFilePath);
%% Begin Parallel Pool
if isempty(gcp) % If there is no current parallel pool
    parpool(NUM_WORKERS) % Create parallel pool on default cluster of size NUM_WORKERS
end
%% Data Processing: Load data using load_data() function
[ data, regressionModel, domainBoundaries, predictionVector, theta, varEps ] = load_data(dataSource, nXGrid, nYGrid, offsetPercentage, verbose);
%% Build hierarchical grid structure using build_structure() function
[ knots, ~, nRegions, outputData, predictionLocations ] = build_structure(NUM_LEVELS_M, NUM_PARTITIONS_J, NUM_KNOTS_r, domainBoundaries, offsetPercentage, verbose, data(:,1:3), predictionVector);
%% Switch cases for calculationType
switch calculationType
    case 'optimize'        
        %% Optimize
        isPredicting = false;
        fun = @(thetaOpt)MRA([thetaOpt(1) thetaOpt(2)], outputData, knots, ...
            NUM_LEVELS_M, NUM_PARTITIONS_J, nRegions, isPredicting, verbose, thetaOpt(3));
        % Dummy values required by optimization routine
        A = []; b = []; Aeq = []; beq = [];
        % fmincon() optimizes over the bounds set
        tic; x = fmincon(fun, initalEstimate, A, b, Aeq, beq, lowerBound, upperBound);
        optRunTime = toc;  % Unsuppress output to print to command window
        % Assign values from optimization to theta and varEps
        theta = [x(1) x(2)];
        varEps = x(3);
        save([resultsFilePath, 'Optimization_Results'], 'theta', 'varEps');
    case 'prediction'
        %% Prediction
        isPredicting = true;
        tic;
        [ sumLogLikelihood, predictions ] = MRA(theta, outputData, knots, ...
            NUM_LEVELS_M, NUM_PARTITIONS_J, nRegions, isPredicting, verbose, varEps, predictionLocations);
        elapsedTimePrediction = toc;  % Unsurpress output to print to command window
        % Reformat data for plotting
        predictions = cell2mat(predictions); predictionLocations = cell2mat(predictionLocations); predictionVariance = predictions(:,2);
        % Add the prediction from the regression
        predRegression = predict(regressionModel, predictionLocations);
        predictionMean = predictions(:,1) + predRegression;
        save([resultsFilePath,'Prediction_Results_', dataSource], 'predictionLocations', 'predictionMean', 'predictionVariance');
        %% Plots
        if displayPlots || savePlots % If plotting
            create_plots(data, predictionLocations, predictionMean, predictionVariance, verbose, displayPlots, savePlots, plotsFilePath)
        end
    case 'likelihood'
        %% Likelihood
        isPredicting = false;
        tic;
        [ sumLogLikelihood ] = MRA(theta, outputData, knots, ...
            NUM_LEVELS_M, NUM_PARTITIONS_J, nRegions, isPredicting, verbose, varEps);  % Unsuppress output to print to command window
        elapsedTimeLikelihood = toc; % Unsuppress output to print to command window
        if verbose % Display the sumLogLikelihood
            sumLogLikelihood
        end
        save([resultsFilePath, 'Likelihood_Results'], 'sumLogLikelihood');
    otherwise
        error('Undefined calculationType. Code is not executed.')
end
if verbose
    disp('MRA execution completed');
end
