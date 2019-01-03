function [ predictionsj ] = spatial_prediction( posteriorPredictionMeanj, posteriorPredictionVariancej, Btildej, RposteriorCholb, KcAb, Kcwb )
%% SPATIAL_PREDICTION.m calculates the posterior predictive distribution
%
%	Input: posteriorMeanj, posteriorVariancej, Btildej, RposteriorCholb, KcAb, Kcwb
%
%	Output: predictionsj
%
%%

% Prediction at finest resolution at the end
NUM_LEVELS_M = size(posteriorPredictionMeanj,2); % helper construct to get M
% Loop through iLevel, making posterior basis function matrices
% Number of basis function matrices is a function of level
for iLevel = NUM_LEVELS_M-1 : -1 : 1    
    KcBtilde = RposteriorCholb{iLevel} \ Btildej{iLevel+2}{iLevel}';    
    Btildej{iLevel+1} = cell(iLevel,1);    
    for jLevel = (iLevel-1) : -1 : 1
        % Equation 13, "posterior basis-function matrices"
        Btildej{iLevel+1}{jLevel} = Btildej{iLevel+2}{jLevel} - KcBtilde' * KcAb{iLevel}{jLevel};
    end
end

for iLevel = 1 : NUM_LEVELS_M-1   
    KcBtildeCurrent = RposteriorCholb{iLevel}\Btildej{iLevel+2}{iLevel}';
    % Equation 12, posterior mean for the weights (eta's)
    posteriorPredictionMeanj(:,iLevel) = KcBtildeCurrent' * Kcwb{iLevel};
    % Equation 12, posterior variance for the weights (eta's)
    posteriorPredictionVariancej(:,iLevel) = sum(KcBtildeCurrent.^2,1);
end

predictionsj= [sum(posteriorPredictionMeanj,2),sum(posteriorPredictionVariancej,2)];

end