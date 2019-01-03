function [ xMin, xMax, yMin, yMax ] = create_partition( xMin0, xMax0, yMin0, yMax0, NUM_PARTITIONS_J )
%% CREATE_PARTITION specifies child partitions by creating a linear space
%   in both the x-direction and the y-direction with J breaks.    
%
%   Partitions are specified by their initial minimal and maximal x and y values
% 
%   Input: xMin0 (double), xMax0 (double), yMin0 (double), yMax0 (double),
%   NUM_PARTITIONS_J (double)
%
%   Output: xMin (double), xMax (double), yMin (double), yMax (double) 
%   minimal and maximal x and y values for j^th region at a particular level 
%%
if ~(NUM_PARTITIONS_J == 2 || NUM_PARTITIONS_J == 4) % If J is not 2 or 4 
    error('myfuns:create_partition: J must be either 2 or 4')
end

switch  NUM_PARTITIONS_J
    % Case for J = 2
    case 2
        % If the X dimension is larger or both are equal
        if (xMax0-xMin0) >= (yMax0-yMin0) 
            tempx = linspace(xMin0,xMax0, NUM_PARTITIONS_J+1);
            xMin = repmat(tempx(1:end-1), NUM_PARTITIONS_J/2, 1); xMin = xMin(:);
            xMax = repmat(tempx(2:end), NUM_PARTITIONS_J/2, 1); xMax = xMax(:);
            
            tempy = linspace(yMin0,yMax0, NUM_PARTITIONS_J/2+1);
            yMin = repmat(tempy(1:end-1),1, NUM_PARTITIONS_J/2+1); yMin = yMin';
            yMax = repmat(tempy(2:end),1, NUM_PARTITIONS_J/2+1); yMax = yMax';
            
        else
            tempx = linspace(xMin0,xMax0, NUM_PARTITIONS_J/2+1);
            xMin = repmat(tempx(1:end-1),1, NUM_PARTITIONS_J/2+1); xMin = xMin(:);
            xMax = repmat(tempx(2:end),1, NUM_PARTITIONS_J/2+1); xMax = xMax(:);
            
            tempy = linspace(yMin0,yMax0, NUM_PARTITIONS_J+1);
            yMin = repmat(tempy(1:end-1),1, NUM_PARTITIONS_J/2); yMin = yMin';
            yMax = repmat(tempy(2:end),1, NUM_PARTITIONS_J/2); yMax = yMax';        
        end
    % Case for J = 4
    case 4
        tempx = linspace(xMin0,xMax0, NUM_PARTITIONS_J/2+1);
        xMin = repmat(tempx(1:end-1), NUM_PARTITIONS_J/2, 1); xMin = xMin(:);
        xMax = repmat(tempx(2:end), NUM_PARTITIONS_J/2, 1); xMax = xMax(:);
        
        tempy = linspace(yMin0,yMax0, NUM_PARTITIONS_J/2+1);
        yMin = repmat(tempy(1:end-1),1, NUM_PARTITIONS_J/2); yMin = yMin';
        yMax = repmat(tempy(2:end),1, NUM_PARTITIONS_J/2); yMax = yMax';        
end
end