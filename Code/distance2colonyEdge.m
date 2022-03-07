function varargout = distance2colonyEdge(colony, line, orientation, varargin)
% Calculate the minimal distance of some points (in 2D space or more)
% towards one curve (sapce with the same dimension). Here coarse graining
% technology is used to accelerate the computation.
%
% Usage: Dis = distance2colonyEdge(colony, line, orientation)
%        distance2colonyEdge(colony, line, orientation)
%        Dis = distance2colonyEdge(..., 'edgeWidth', 10)
%        distance2colonyEdge(..., 'edgeWidth', 10)
%
% COLONY should be a binary matrix, the logical value 1 represent the
% pixels of biofilm.
%
% LINE should be a 2*2 matrix, specifying the boundary line used to remove
% the 'pseudo edge' contacting the chamber wall. The format is [x1, y1; x2,
% y2]
%
% ORIENTATION should be a string, either 'above' or 'below', specifying the
% position of 'pseudo edge' relative to the line.
%
% The optional input EDGEWIDTH determines the definition of the edge
% region.region. The values should be in units of pixels.
%
% Output
% The output DISTANCE is a n*1 array, presenting the minimal distance of n
% points towards the curve.
%
% Written by Tianmin Wang @ Tsinghua University
% Version 0.0.1. Created on Aug 08, 2019. Last modified on Aug 08, 2019.

argin = inputParser;
addOptional(argin,'edgeWidth',3)
parse(argin,varargin{:})
edgeWidth = argin.Results.edgeWidth;

if ~islogical(colony)
    error('COLONY must be a logical array.')
end

[Ny, Nx] = size(line);
if Ny ~= 2 || Nx ~= 2
    error('LINE should be a 2*2 array!')
end

if ~isempty(find(contains({'above', 'below'}, orientation), 1))
    aboveFlag = strcmp('above', orientation);
else
    error('ORIENTATION should be "above" or "below"!')
end

% calculate colony edge
colonyE = colonyEdge(colony, edgeWidth);

% remove edge between biofilm and chamber boundary (mode 2)
% 1. calculate the line defined by LINE
coefficients = polyfit([line(1,1), line(2,1)], ...
                       [line(1,2), line(2,2)], 1);
slope = coefficients (1);
intercept = coefficients (2);

% 2. convert mode 2 edge to non-edge
colonyE_v2 = colonyE;
[Ny, Nx] = size(colonyE_v2);
for y = 1:Ny
    for x = 1:Nx
        coorFlag = (y - (slope*x + intercept)) > 0; % whether the pixel is above or below the line
        if coorFlag
            if ~aboveFlag
                colonyE_v2(y, x) = 0;
            end
        else
            if aboveFlag
                colonyE_v2(y, x) = 0;
            end
        end
    end
end

% calculate the distance to edge
dis2edge = bwdist(colonyE_v2);

% output
if nargout == 0
    % present the histogram of each channel
    contour(dis2edge);
    colorbar;
else
    varargout{1} = dis2edge;
    varargout{2} = colonyE_v2;
end

