function [ M_integer, nRegionsFinest, M_estimate, nRegions, totalRegions ] = find_num_levels_suggested( NUM_DATA_POINTS_n, NUM_KNOTS_r, NUM_PARTITIONS_J )
%% FIND_NUM_LEVELS_SUGGESTED is a stand-alone function that estimates number of levels
%   Finds the number of levels for a given number of observations, number
%   of knots per region (NUM_KNOTS_r) and number of partitions per region (NUM_PARTITIONS_J).
%   The idea is to work backwards from the finest level by making the average
%   number of observations per region similar to the number of knots.
%
%	Input: NUM_DATA_POINTS_n (double), NUM_KNOTS_r (double),
%	NUM_PARTITIONS_J (double)
%
%	Output: M_integer (double), nRegionsFinest (double), M_estimate (double),
%	nRegions (vector), totalRegions (double)
%%
% Find number of regions required at finest level
nRegionsFinest = NUM_DATA_POINTS_n/NUM_KNOTS_r;

% Solve for M which is m at the finest level using the J^M is the number of
% regions at the finest level

M_estimate = log(nRegionsFinest)/log(NUM_PARTITIONS_J);
M_integer = int64(ceil(M_estimate));

% Comment: Matlab doesn't have a direct command to solve a logarithm for
% any base, so we use log(x)/log(b) in place of log_base_b(x)

mLevels = 0:M_integer; % Vector of levels
nRegions = NUM_PARTITIONS_J.^mLevels; % Regions (partitions) by level
totalRegions = int64(sum(nRegions)); % Total number of regions over all levels
end