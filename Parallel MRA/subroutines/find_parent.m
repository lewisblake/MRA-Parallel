function [ levelParent, tileParent, indexParent ] = find_parent( index, nRegions, NUM_PARTITIONS_J )
%% FIND_PARENT finds descriptors of parent
%   This function finds the level, tile number and continuous index of the
%   parent of the region for the given index
%
%   Input: index, nRegions, NUM_PARTITIONS_J
%
%   Output: levelParent, tileParent, indexParent
%%
switch  index
    case 1
        % Special case for zeroth region
        tileParent = NaN;
        levelParent = NaN;
        indexParent = NaN;
    otherwise
        cumulativeRegions = cumsum(nRegions); % vector of the cumulative number of regions up until each level
        indexSmaller = find(cumulativeRegions < index, 1, 'last'); % cumulativeRegions will have largest index at coarser resolution
        level = indexSmaller + 1; % Present level
        tile = index - cumulativeRegions(indexSmaller);
        % Figure out ID of parent
        tileParent = ceil(tile/NUM_PARTITIONS_J);
        levelParent = level-1;
        indexParent = sum(nRegions(1 : level-2)) + tileParent;
end
end