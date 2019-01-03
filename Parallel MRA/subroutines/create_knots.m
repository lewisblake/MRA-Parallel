function [ X, Y ] = create_knots( xMin, xMax, nX, yMin, yMax, nY, offsetPercentage )
%% CREATE_KNOTS creates the knots for each partition at each level 
%   (Determines knot locations for each partition)
%
%	Input: xMin (double), xMax (double), nX (double), yMin (double), yMax
%	(double), nY (double), offsetPercentage (double)
%	
%	xmin/xmax are the mimimum and maximum x-values, similarly for ymin/ymax
%	nx is the number of knots in the x-direction, similarly for ny
%	offsetPercentage is how much (as a percentage), the spatial domain used throughout calculations
%	is smaller than the maximum domain as determined by xmin/xmax and ymin/ymax
%
%	Output: [X, Y], which is a meshgrid of the knots in X-Y space (matrix)
%%
offsetX = (xMax-xMin)*offsetPercentage/100;
offsetY = (yMax-yMin)*offsetPercentage/100;
[X,Y] = meshgrid(linspace(xMin + offsetX, xMax - offsetX, nX), linspace(yMin + offsetY, yMax - offsetY, nY));
end