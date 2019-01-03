function [ sumLogLikelihood, predictions ] = MRA( theta, data, knots, ...
    NUM_LEVELS_M, NUM_PARTITIONS_J, nRegions, isPredicting, verbose,  varargin )
%% MRA.m is the primary Multi-resolution approximation function
%
%  First, this function pre-allocates the memory needed for objects used in
%  subsequent calculations. Then, the prior distribution is created by
%  looping from coarsest resolution to finest resolution, where at each
%  iLevel, calculations are completed in parallel with use of parfor. After
%  the prior has been created, the algorithm works from second-finest
%  resolution to coarsest resolution computing the posterior inference.
%  Again, the algorithm loops between iLevel's in serial where at each
%  iLevel, computations are computed in parallel with use of parfor.
%  Finally, if in prediction mode, spatial prediction takes with
%  with the use of the spatial_prediction() function in parallel.
%
%  Input: theta (vector), data (cell), knots (cell), NUM_LEVELS_M (double),
%  NUM_PARTITIONS_J (double), nRegions (vector), isPredicting (boolean),
%  verbose (boolean), varargin (varargin can be varEps(double) and predictionLocationsj (matrix))
%
%  Output: sumLogLikelihood (double), predictions (matrix)
%%
% Check number of optional input arguments
numVarArgsIn = length(varargin);
if numVarArgsIn > 2
    error('myfuns:create_prior:TooManyInputs', ...
        'requires at most 2 optional inputs');
end
% Set defaults for optional inputs
optionalArguments = {0 NaN};
% overwrite the ones specified in varargin.
optionalArguments(1 : numVarArgsIn) = varargin;
[ varEps, predictionLocationsj ] = optionalArguments{:};
% Calculate key quantities
totalRegions = sum(nRegions);
cumulativeRegions = cumsum(nRegions);
% Pre-allocate quantities
RpriorChol = cell(totalRegions,1);
KcB = cell(totalRegions,1);
AtildePrevious = cell(totalRegions,1);
wtildePrevious = cell(totalRegions,1);
logLikelihood = NaN(totalRegions,1);
if isPredicting % If predicting, pre-allocate space for necessary quantities
    posteriorPredictionMean = cell(nRegions(NUM_LEVELS_M),1); % Only need to store values at finest resolution
    posteriorPredictionVariance = cell(nRegions(NUM_LEVELS_M),1); % Only need to store values at finest resolution
    Btilde = cell(nRegions(NUM_LEVELS_M),1); % Only need to store values at finest resolution
    predictions = cell(nRegions(NUM_LEVELS_M),1); % Only need to store values at finest resolution
    RposteriorChol = cell(totalRegions,1);
    KcholA = cell(totalRegions,1);
    Kcholw = cell(totalRegions,1);
else % Workaround to pass correct object to create_prior
    predictionLocationsj = num2cell(nan(totalRegions,1));
end

%% Prior distribution
% Preallocate cell containers for create_prior()
KcBContainer = cell(NUM_LEVELS_M-1,nRegions(NUM_LEVELS_M));
RpriorCholContainer = cell(NUM_LEVELS_M-1,nRegions(NUM_LEVELS_M));
knotsContainer = cell(NUM_LEVELS_M,nRegions(NUM_LEVELS_M)); % Rows include knots at this iLevel
% Loop from coarsest to finest level
for iLevel = 1 : NUM_LEVELS_M
    % Make vector of indices for this level
    dummyIndexVec = (cumulativeRegions(iLevel) - nRegions(iLevel) + 1) : cumulativeRegions(iLevel);
    % Find "number of indices for this iLevel" = nPartitionsThisLevel
    nRegionsAtThisLevel = length(dummyIndexVec); %
    % Determine continuous the first continous index for the level (iLevel)
    if iLevel == 1
        indexBeginningThisLevel = 0;
    else
        indexBeginningThisLevel = cumulativeRegions(iLevel) - nRegions(iLevel);
    end
    % Create Containers for this iLevel using a dummy index for regions at
    % this level
    for jDummyIndex = 1:nRegionsAtThisLevel
        % Continous index for this region
        thisDummyIndex = dummyIndexVec(jDummyIndex);
        % For each region's index in this iLevel, find it's ancestry
        thisIndexAncestry = find_ancestry(thisDummyIndex, nRegions, NUM_PARTITIONS_J);       
        % Populate Containers for this iLevel
        KcBContainer(1:iLevel-1, jDummyIndex) = KcB(thisIndexAncestry, 1);
        RpriorCholContainer(1:iLevel-1, jDummyIndex) = RpriorChol(thisIndexAncestry,1);
        knotsContainer(1:iLevel, jDummyIndex) = knots([thisIndexAncestry;thisDummyIndex],1);
    end
    % Create pointer variables to section of container needed for iLevel
    thisLevelKcBContainer = KcBContainer(1:iLevel-1,1:nRegionsAtThisLevel);
    thisLevelRpriorCholContainer = RpriorCholContainer(1:iLevel-1,1:nRegionsAtThisLevel);
    thisLevelKnotsContainer = knotsContainer(1:iLevel, 1:nRegionsAtThisLevel);   
    % Parallel loop through all indices for a given iLevel
    parfor jIndex = 1:nRegionsAtThisLevel
        % Calculate the prior quantities for this region using create_prior()
        [thisRpriorChol, thisKcholBchol, thisAtj, thiswtj, thisRetLikPred] = create_prior(theta, ...
            NUM_LEVELS_M, thisLevelKnotsContainer(:,jIndex), thisLevelRpriorCholContainer(:,jIndex), ...
            thisLevelKcBContainer(:,jIndex), data{indexBeginningThisLevel+jIndex,1}, ...
            varEps, predictionLocationsj{indexBeginningThisLevel+jIndex,1});

        % Collect prior quantities for this region
        RpriorChol{indexBeginningThisLevel+jIndex} = thisRpriorChol;
        KcB{indexBeginningThisLevel+jIndex} = thisKcholBchol;       
        % If at finest iLevel, store quantities for posterior inference
        if iLevel == NUM_LEVELS_M
            AtildePrevious{indexBeginningThisLevel+jIndex} = thisAtj;
            wtildePrevious{indexBeginningThisLevel+jIndex} = thiswtj;
            if isPredicting  % If predicting
                posteriorPredictionMean{jIndex} = thisRetLikPred{1};
                posteriorPredictionVariance{jIndex} = thisRetLikPred{2};
                Btilde{jIndex} = thisRetLikPred{3};
            else % If not predicting
                logLikelihood(indexBeginningThisLevel+jIndex,1) = thisRetLikPred;
            end
        end
    end
    if verbose
        disp(['Prior Level ',num2str(iLevel), ' completed'])
    end
end
% Uncomment below if needed to save memory. Setting to empty vectors forces freeing memory as opposed to using clear.
%KcB = []; KcBContainer = []; RpriorCholContainer = []; knotsContainer = []; knots = [];

%% Posterior inference
% Pre-allocate posterior quantities
AtildeCurrent = cell(totalRegions,1);
wtildeCurrent = cell(totalRegions,1);

% Loop iLevel from second finest to coarsest resolution
for iLevel = NUM_LEVELS_M-1 : -1 : 1
    if verbose
        disp(['Posterior Level ',num2str(iLevel),' starting']);
    end
    % For all levels except coarsest one, set begining index
    if iLevel == 1
        indexBeginningThisLevel = 0;
    else
        indexBeginningThisLevel = cumulativeRegions(iLevel) - nRegions(iLevel);
    end
    % From the entire cell wtildePrevious, find those cell entries that
    % correspond to the children of the current level. Cells are ordered.
    wtildeFinerLevel = wtildePrevious(NUM_PARTITIONS_J^iLevel: NUM_PARTITIONS_J^(iLevel + 1) - 1);
    AtildeFinerLevel = AtildePrevious(NUM_PARTITIONS_J^iLevel: NUM_PARTITIONS_J^(iLevel + 1) - 1);
    % Reshape wtildeFinerLevel & AtileFinerLevel to have each subregion's
    % values in columns. Each tile had J children at finer level.
    % reshape() is effectively zero-cost function
    wtildeFinerLevel = reshape(wtildeFinerLevel, NUM_PARTITIONS_J, []);
    AtildeFinerLevel = reshape(AtildeFinerLevel, NUM_PARTITIONS_J, []);    
    % Set number of partitions for this level
    nRegionsAtThisLevel = length((cumulativeRegions(iLevel) - nRegions(iLevel) + 1) : cumulativeRegions(iLevel));
    % Parallel loop through all regions for the level
    parfor jRegion = 1:nRegionsAtThisLevel
        % Set prior covariance for this tile
        RpriorCholj = RpriorChol{indexBeginningThisLevel+jRegion};
        % Calculate posterior inference using posterior_inference()
        % Column indexing into wtildeFinerLevel, etc. due to reshape function
        [ wtildeCurrentj, AtildeCurrentj, logLikelihoodj, ...
            RposteriorCholj, Kcholwj, KcholAj ] = posterior_inference(iLevel, RpriorCholj,...
            wtildeFinerLevel(:, jRegion), AtildeFinerLevel(:, jRegion));
        
        wtildeCurrent{indexBeginningThisLevel+jRegion} = wtildeCurrentj;
        AtildeCurrent{indexBeginningThisLevel+jRegion} = AtildeCurrentj;
        
        if isPredicting
            RposteriorChol{indexBeginningThisLevel+jRegion} = RposteriorCholj;
            Kcholw{indexBeginningThisLevel+jRegion} = Kcholwj;
            KcholA{indexBeginningThisLevel+jRegion} = KcholAj;
        else
            logLikelihood(indexBeginningThisLevel+jRegion,1) = logLikelihoodj;
        end
    end
    wtildePrevious = wtildeCurrent;
    AtildePrevious = AtildeCurrent;
end
% Uncomment below if needed to save memory. Setting to empty vectors forces freeing memory as opposed to using clear.
% wtildePrevious = []; AtildePrevious = []; RpriorChol = []; % To save memory
sumLogLikelihood = sum(logLikelihood);

%% Spatial prediction
if isPredicting
    if (NUM_LEVELS_M > 0)
        % Vector of indices for finest level, NUM_LEVELS_M
        dummyIndexVec = (cumulativeRegions(NUM_LEVELS_M) - nRegions(NUM_LEVELS_M) + 1) : cumulativeRegions(NUM_LEVELS_M);
        % Number of partitions at finest level, NUM_LEVELS_M
        nRegionsAtThisLevel = length(dummyIndexVec);
        % Begining of the continous index for finest level, NUM_LEVELS_M
        indexBeginningThisLevel = cumulativeRegions(NUM_LEVELS_M) - nRegions(NUM_LEVELS_M);
        % Containers for finest level, NUM_LEVELS_M
        RposteriorCholContainer = cell(NUM_LEVELS_M-1, nRegionsAtThisLevel);
        KcAContainer = cell(NUM_LEVELS_M-1, nRegionsAtThisLevel);
        KcwContainer = cell(NUM_LEVELS_M-1, nRegionsAtThisLevel);
        % Fill containers for this level
        for jDummyIndex = 1:nRegionsAtThisLevel
            dummyIndex = dummyIndexVec(jDummyIndex);
            thisIndexAncestry = find_ancestry(dummyIndex, nRegions, NUM_PARTITIONS_J);
            RposteriorCholContainer(:,jDummyIndex) = RposteriorChol(thisIndexAncestry,1);
            KcAContainer(:,jDummyIndex) = KcholA(thisIndexAncestry,1);
            KcwContainer(:,jDummyIndex) = Kcholw(thisIndexAncestry,1);
        end
        % Parallel loop over partitions at finest level making predicitions
        parfor jRegion = 1:nRegionsAtThisLevel
            predictions{indexBeginningThisLevel+jRegion,1} = spatial_prediction(posteriorPredictionMean{jRegion,1},...
                posteriorPredictionVariance{jRegion,1}, Btilde{jRegion,1}, ...
                RposteriorCholContainer(:,jRegion), KcAContainer(:,jRegion), KcwContainer(:, jRegion));
        end        
    else
        dummyIndex = 1; % Dummy counter to index into predictions, etc. correctly
        for index = (cumulativeRegions(NUM_LEVELS_M) - nRegions(NUM_LEVELS_M) + 1) : cumulativeRegions(NUM_LEVELS_M)
            predictions{dummyIndex,:} = [posteriorPredictionMean(dummyIndex,1), posteriorPredictionVariance(dummyIndex,1)];
            dummyIndex = dummyIndex + 1;
        end
    end
end