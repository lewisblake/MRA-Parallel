function [ knots, partitions, nRegions, outputData, predictionLocations ] = build_structure( NUM_LEVELS_M, ...
    NUM_PARTITIONS_J, NUM_KNOTS_r, domainBoundaries, offsetPercentage, verbose, varargin )
%% BUILD_STRUCTURE builds the nested partitioning structure as a function of levels
%	Input: NUM_LEVELS_M (double), NUM_PARTITIONS_J (double), NUM_KNOTS_r (double), 
%   domainBoundaries (vector), offsetPercentage (double), varargin (varargin allows inputting of data (matrix) and the predictionVector (matrix))
%	
%	NUM_LEVELS_M = Total number of levels (levels are indexed by level)
%	NUM_PARITIONS_J = Number of partitions for each region for a given level
%	NUM_KNOTS_r = Number of knots in each partition
% 	domainBoundaries is a vector containing the minimal and maximal x and y values
%
%	Output:
% 	build_structure returns the knots (cell), partitions (cell), nRegions (vector), 
%	outputData (cell), and the predictionLocations (matrix)
%%
% Check number of optional input arguments does not exceed two
numVarArgs = length(varargin);
if numVarArgs > 2
    error('myfuns:build_structure:TooManyInputs', ...
        'requires at most 2 optional inputs');
end

% Error handling on input arguments
if ~isinf(NUM_LEVELS_M) && floor(NUM_LEVELS_M) ~= NUM_LEVELS_M
    error('myfuns:build_structure: NUM_LEVELS_M must be a positive integer.')
elseif ~isinf(NUM_PARTITIONS_J) && floor(NUM_PARTITIONS_J) ~= NUM_PARTITIONS_J
    error('myfuns:build_structure: NUM_PARTITIONS_J must be either 2 or 4.')
elseif ~isinf(NUM_KNOTS_r) && floor(NUM_KNOTS_r) ~= NUM_KNOTS_r
    error('myfuns:build_structure: NUM_KNOTS_r must be a positive integer.')
elseif offsetPercentage < 0 || offsetPercentage >= 1
    error('myfuns:build_structure: offsetPercentage must be between 0 and 1.')
end

% Optional argument that can be passed is data
optArgs(1:numVarArgs) = varargin;
[ data, predictionVector ] = optArgs{:};

%% Calculate quantities of interest
mLevels = 0:NUM_LEVELS_M-1; % Create a vector of levels
nRegions = NUM_PARTITIONS_J.^mLevels; % Vector of regions (partitions) at each level
totalRegions = sum(nRegions); % Calculate total number of regions

% Calculate number of knots in each direction
if isinteger(sqrt(NUM_KNOTS_r))  % Assign knots.
nKnotsX0=sqrt(NUM_KNOTS_r); nKnotsX=sqrt(NUM_KNOTS_r); % Number of knots in x-direction
nKnotsY0=sqrt(NUM_KNOTS_r); nKnotsY=sqrt(NUM_KNOTS_r); % Number of knots in y-direction
else
nKnotsX0=ceil(sqrt(NUM_KNOTS_r)); nKnotsX=ceil(sqrt(NUM_KNOTS_r));
nKnotsY0=NUM_KNOTS_r/nKnotsX0; nKnotsY=NUM_KNOTS_r/nKnotsX;  
end
%% Create knots for partitions
% Pre-allocate
knots = cell(totalRegions,1);
outputData = cell(totalRegions,1);
partitions = cell(totalRegions,1);

% Construct zeroth level
xMin0 = domainBoundaries(1); xMax0 = domainBoundaries(2);
yMin0 = domainBoundaries(3); yMax0 = domainBoundaries(4);
% Edge buffer added to xMax0 and yMax0 to
% include all observations on the boundary at the zeroth level
xMax0 = xMax0 + (offsetPercentage/2)*(xMax0 - xMin0);
yMax0 = yMax0 + (offsetPercentage/2)*(yMax0 - yMin0);

[ knotsX, knotsY ] = create_knots(xMin0, xMax0, nKnotsX0, yMin0, yMax0, nKnotsY0, offsetPercentage);
knots{1,1}=[knotsX(:), knotsY(:)];
[ xMin, xMax, yMin, yMax ] = create_partition(xMin0, xMax0, yMin0, yMax0, NUM_PARTITIONS_J);
partitions{1,1}=[ xMin, xMax, yMin, yMax ];

% Set finest knot level
if numVarArgs == 2
    finestKnotLevel = NUM_LEVELS_M-1;
else
    finestKnotLevel = NUM_LEVELS_M;
end

% Loop through each level
% Begin at level 2 since zeroth level is indexed as first level
for iLevel = 2:finestKnotLevel 
    if verbose % Display progress indicators
        disp(['Building Level ',num2str(iLevel),' starting']);
    end
    % At each level, create tiles for partitioning into J subregions
    % Loop through tiles at each resolution
    for jTile = 1:NUM_PARTITIONS_J:nRegions(iLevel)
        % Find continuous index ID (i) of parent
        [ index ] = find_index( iLevel, jTile, nRegions );
        [ ~, ~, indexParent ] = find_parent(index, nRegions, NUM_PARTITIONS_J);
        % Get partition coordinates of parent
        xMin = partitions{indexParent,1}(:,1);
        xMax = partitions{indexParent,1}(:,2);
        yMin = partitions{indexParent,1}(:,3);
        yMax = partitions{indexParent,1}(:,4);
        
        for kPartition = 1:NUM_PARTITIONS_J
            indexCurrent = sum(nRegions(1:iLevel-1)) + jTile + kPartition - 1;            
            [knotsX,knotsY] = create_knots(xMin(kPartition), xMax(kPartition), nKnotsX, yMin(kPartition), yMax(kPartition), nKnotsY, offsetPercentage);           
            [ xMinTemp, xMaxTemp, yMinTemp, yMaxTemp ] = create_partition(xMin(kPartition), xMax(kPartition), yMin(kPartition), yMax(kPartition), NUM_PARTITIONS_J);           
            knots{indexCurrent,1} = [knotsX(:),knotsY(:)];
            partitions{indexCurrent,1} = [ xMinTemp, xMaxTemp, yMinTemp, yMaxTemp ];            
        end
    end
end

% Special construct to find knots at the finest resolution level
if numVarArgs == 2  % The last level is built using data as knots
    iLevel = NUM_LEVELS_M; % Set the level as the finest level
    if verbose % Display progress indicators
        disp(['Building finest resolution Level ',num2str(iLevel),' starting']);
    end
    nTilesFinestLevel = nRegions(end)/NUM_PARTITIONS_J;    
    % Loop through tiles at each resolution
    for jTile = 1:NUM_PARTITIONS_J : nRegions(iLevel)       
        % Figure out parent index
        [ index ] = find_index(iLevel, jTile, nRegions);
        [ ~, ~, indexParent ] = find_parent(index, nRegions, NUM_PARTITIONS_J);
        % Get partition coordinates of parent
        xMin = partitions{indexParent,1}(:,1);
        xMax = partitions{indexParent,1}(:,2);
        yMin = partitions{indexParent,1}(:,3);
        yMax = partitions{indexParent,1}(:,4);
        % Loop through partitions at each tile
        for kPartition = 1:NUM_PARTITIONS_J
            indexCurrent = sum(nRegions(1:iLevel-1))+jTile+kPartition-1;            
            ind = find(data(:,1) >= xMin(kPartition) & data(:,1) < xMax(kPartition) & data(:,2) >= yMin(kPartition) & data(:,2) < yMax(kPartition));            
            knotsX = data(ind,1); knotsY = data(ind,2);            
            knots{indexCurrent,1} = [knotsX(:),knotsY(:)];            
            outputData{indexCurrent,1} = data(ind,3); % This is to only pass the data, not the location to MRA            
            data(ind,:) = []; % Eliminate the data that has already been assigned to a region, speeds up subsequent searching  
            %% Vinay's addition to partition the prediction locations
            if ~isnan(predictionVector)  % If predicting
            predictionIndex = find(predictionVector(:,1) >= xMin(kPartition) & predictionVector(:,1) < xMax(kPartition) & predictionVector(:,2) >= yMin(kPartition) & predictionVector(:,2) < yMax(kPartition));
            predictionLocations{indexCurrent,1} = predictionVector(predictionIndex,:);
            else
                predictionLocations = NaN;
            end            
        end
        if verbose % Display progress indicators
            disp(['Tile ',num2str(floor(jTile/NUM_PARTITIONS_J)+1),' of ', num2str(nTilesFinestLevel), ' completed']);
        end
    end
else
end
end
