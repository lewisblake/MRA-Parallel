function [ index ] = find_index( level, tileNum, nRegions )
%% FIND_INDEX finds continuous index
%   This function finds the continuous index given the level and
%   tile as inputs
%   
%   index is "continuous" in the sense that it takes all integer values from
%   1 to totalRegions
%
%	Input: level, tileNum, nRegions
%
%	Output: index, continous index 
%%
index = sum(nRegions(1 : level-1)) + tileNum;
end