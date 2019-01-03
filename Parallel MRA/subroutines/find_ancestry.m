function [ indexAncestry ] = find_ancestry( index, nRegions, NUM_PARTITIONS_J )
%% FIND_ANCESTRY finds ancestors of a region
%   This function finds the level, tile number and continuous index of the
%   ancestry of the region with given index
%
%   Input: index, nRegions, NUM_PARTITIONS_J
%
%   Output: indexAncestry, indexes for parental hierarchy (ancestry)
%%
switch  index
    case 1
        % Special case for zeroth region
        indexAncestry = [];
    
    otherwise
        cummulativeRegions = cumsum(nRegions);
        indexSmaller = find(cummulativeRegions < index, 1, 'last');
        % pre-allocate vector
        indexAncestry = zeros(indexSmaller, 1, 'int64'); % number of ancestry members is one smaller than level         
         % Fill ancestry by looping through parent's parents
         for k = 1 : indexSmaller
            [ ~,~,indexParent ] = find_parent( index, nRegions,NUM_PARTITIONS_J );           
            indexAncestry(k,1) = indexParent;           
            index = indexParent;
         end
         indexAncestry = flip(indexAncestry); % Order from coarsest to finest
end
end