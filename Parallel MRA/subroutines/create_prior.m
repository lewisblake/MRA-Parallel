function [ RpriorCholj, KcBc, Atj, wtj, retLikelihoodPred ] = create_prior( theta, ...
    NUM_LEVELS_M, knotsb, RpriorCholb, KcBb, dataj, varargin )
%% CREATE_PRIOR creates prior values
%   for current region and ancestry; this function contains optional input
%   and output arguments depending on whether the level is the last level
%   and/or if predictions are made, optional outputs are Atj,wtj,Btildej ("A tilde j", etc.)
%
%   Input: theta, M, knotsb, RpriorCholb, KcBb, dataj, varargin
%
%   Output: RpriorCholj, KcBc, Atj, wtj, retlikpred
%%

% Check number of optional input arguments
numvarargsin = length(varargin);
if numvarargsin > 2
    error('myfuns:create_prior:TooManyInputs', ...
        'requires at most 2 optional inputs');
end

% Set defaults for optional inputs
optionalArguments = {0 NaN};
% Overwrite the ones specified in varargin.
optionalArguments(1 : numvarargsin) = varargin;
[ varEps, predictionLocationsj ] = optionalArguments{:};


% Preliminaries
currentLevel = length(knotsb); % Find current level
isCurrentLevelLessThanLastLevel = (currentLevel < NUM_LEVELS_M); % Indicator if current level is less than last level

% Create prior quantities
RpriorCholj = [];
KcBc = cell(currentLevel-1, 1);
V = cell(currentLevel, 1);

for iLevel = 1 : currentLevel
    % For each iLevel, compute the covariance matrix between knots at
    % iLevel and knots at currentLevel
    V{iLevel,1} = evaluate_covariance(knotsb{currentLevel,1}, knotsb{iLevel,1}, theta);   
    for k = 1 : (iLevel-1)       
        if iLevel < currentLevel % For every iLevel before the currentLevel            
            % Equation 6, denoted using V instead of W
            V{iLevel} = V{iLevel} - KcBc{k}'*KcBb{iLevel,1}{k,1};            
        else % If iLevel == currentLevel           
            V{iLevel} = V{iLevel} - KcBc{k}'*KcBc{k};
        end
    end
    
    if iLevel < currentLevel % For every iLevel before the currentLevel        
        KcBc{iLevel} = RpriorCholb{iLevel} \ V{iLevel}';      
    else % If iLevel == currentLevel
        % Add diagonal matrix of varEps
        V{iLevel} = V{iLevel} + diag(linspace(varEps,varEps,size(V{iLevel},1)));
        % Compute the Cholesky decompositon for R_prior
        RpriorCholj = chol(V{iLevel}, 'lower');        
    end    
end

% Deal with last level, case L = M, separately
if ~isCurrentLevelLessThanLastLevel % Check if region is at lowest level
    % Begin inference at lowest level
    % Pre-compute solves    
    Sicy = RpriorCholj \ dataj;
    SicB = cell(currentLevel-1, 1);    
    
    for iLevel = 1 : (currentLevel-1)      
        SicB{iLevel} = RpriorCholj \ V{iLevel};       
    end   
    % Inference quantities
    wtj = cell(currentLevel-1, 1);
    Atj = cell(currentLevel-1, currentLevel-1);
    
    for iLevel = 1 : (currentLevel-1)      
        wtj{iLevel} = SicB{iLevel}'*Sicy;       
        for k = iLevel : (currentLevel-1)
            Atj{iLevel,k} = SicB{iLevel}'*SicB{k};
        end
    end
    
    % Check if predicting
    if isnan(predictionLocationsj)  % If NOT predicting
        logLikelihoodj = 2 * sum(log(diag(RpriorCholj))) + Sicy'*Sicy;
        retLikelihoodPred = logLikelihoodj;
        
    else % If predicting
        RpriorChol = [ RpriorCholb; {RpriorCholj} ];
        KcB = [KcBb; {KcBc}];        
        % Calculate Bp and currentLevel
        KcBp = cell(currentLevel, 1);
        Vp = cell(currentLevel, 1);
        
        for iLevel = 1 : currentLevel
            Vp{iLevel} = evaluate_covariance(predictionLocationsj, knotsb{iLevel}, theta);
            for k = 1 : (iLevel-1)
                Vp{iLevel} = Vp{iLevel} - KcBp{k}'*KcB{iLevel,1}{k,1};
            end
            KcBp{iLevel} = RpriorChol{iLevel} \ Vp{iLevel}';
        end
        
        Vpp = evaluate_covariance(predictionLocationsj, predictionLocationsj, theta); % Covariance matrix of prediction locations
        
        for iLevel = 1 : (currentLevel-1)
            Vpp = Vpp - KcBp{iLevel}'*KcBp{iLevel};
        end
        
        % Initialize prediction inference
        posteriorMeanj = NaN(size(predictionLocationsj,1),currentLevel); % prediction mean matrix for all levels
        posteriorVariancej = NaN(size(predictionLocationsj,1),currentLevel); % prediction variance matrix for all levels
        posteriorMeanj(:,currentLevel) = KcBp{currentLevel}'*Sicy;
        posteriorVariancej(:,currentLevel) = diag(Vpp - KcBp{currentLevel}'*KcBp{currentLevel});
        Btildej = cell(currentLevel,1);
        
        for k = 1 : (currentLevel-1)
            Btildej{currentLevel+1}{k} = Vp{k} - KcBp{currentLevel}'*SicB{k};
        end
        retLikelihoodPred = {posteriorMeanj, posteriorVariancej, Btildej};
    end
else  % If NOT lowest level
    wtj = NaN; Atj = NaN; retLikelihoodPred = NaN; % Assign NaNs to outputs
end
end