function [ wtildeCurrentj, AtildeCurrentj, logLikelihoodj, RposteriorCholj, Kcwj, KcAj ] = posterior_inference( thisLevel, RpriorCholj, wtildeChildren, AtildeChildren )
%% POSTERIOR_INFERENCE: This function calculates posterior inference
% From second finest to coarsest resolution
%
%   Input: thisLevel, Rpriorcholj, wtildechildren, Atildechildren
%
%   Output: wtildeCurrentj, AtildeCurrentj, logLikelihoodj, RposteriorCholj, Kcwj, KcAj
%    
%%

w = cell(thisLevel,1);
A = cell(thisLevel,thisLevel); % This is a square cell structure to hold tiles

for iLevel = 1 : thisLevel 
    temp = []; 
    for c = 1:length(wtildeChildren)
        temp = [temp wtildeChildren{c}{iLevel}];
    end
    w{iLevel} = sum(temp,2); % Calculate w at this level. Equation 9
    clear temp; % To save memory
    
    for k = iLevel:thisLevel
        temp = NaN(size(AtildeChildren{1, 1}{1, 1})); % helper construct
        temp(:,:,length(AtildeChildren)) = NaN(size(AtildeChildren{1, 1}{1, 1})); % helper construct
        for c = 1:length(AtildeChildren)
            temp(:,:,c) = AtildeChildren{c}{iLevel,k};
        end
        A{iLevel,k} = sum(temp,3); % Calculate A. Equation 9
        clear temp; % To save memory
        
    end
end

% Calculate Cholesky of K.inv. Equation 8
Rposterior = RpriorCholj * RpriorCholj' + A{thisLevel,thisLevel};
RposteriorCholj = chol(Rposterior,'lower'); % Compute the Cholesky decomposition

% Pre-compute the solves required later
Kcwj = RposteriorCholj\w{thisLevel};
KcAj = cell(thisLevel-1,1);

for iLevel = 1 : thisLevel-1
    KcAj{iLevel} = RposteriorCholj\A{iLevel,thisLevel}';
end

% Calculate w.tilde and A.tilde

if thisLevel == 1 % If at the first level
    wtildeCurrentj = NaN;
    AtildeCurrentj = NaN;
else
    wtildeCurrentj = cell(thisLevel-1,1); % Each cell holds vectors
    AtildeCurrentj = cell(thisLevel-1,thisLevel); % This is a square cell structure to hold tiles
    for iLevel = 1 : thisLevel-1
        wtildeCurrentj{iLevel} = w{iLevel}-KcAj{iLevel}'*Kcwj;
        for k = iLevel : thisLevel-1
            AtildeCurrentj{iLevel,k} = A{iLevel,k}-KcAj{iLevel}'*KcAj{k};    
        end
    end   
end

logLikelihoodj = 2 * sum(log(diag(RposteriorCholj))) - 2 * sum(log(diag(RpriorCholj))) - Kcwj' * Kcwj;

end