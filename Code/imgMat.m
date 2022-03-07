function varargout = imgMat(M, distance2Edge, channelsOfInterest, colony, varargin)

% This function is to process biofilm img data to a n*m array, ready for
% PCA analysis and clustering;
% n: bacterial cells; m: number of channels
% Usage: [myDataNor, myData] = imgMat_new(M, distance2Edge, channelsOfInterest, colony)
%        imgMat(M, distance2Edge, channelsOfInterest, colony)
%        imgMat(..., 'outlierThr',[0.001,0.001,0.001])
%
% Input: 
% M must be a cell that contain data from all fluorescent channels. Note
% that the order of channels here should be consistent with
% channelsOfInterest. Thus, length(M) = length(channelsOfInterest);
%
% DISTANCE2EDGE is a matrix with the same dimension of M; specifying the
% distance of each colony pixel towards colony edge.
%
% The input channelsOfInterest is a cell array, specifying the channels of
% interest. Accepted channels include 'CFP', 'GFP', 'mCherry' and 'APC'
%
% Colony should be a binary matrix, the logical value 1 represent the
% pixels of biofilm.
%
% Center should be a array with 2 elements represent x and y coordiante.
% Using iViewer(M, 'roi', @impoint) can get the biofilm center coordinate.
%
% The optional input outlierThr determines largest proportion of data in
% one channel, that is regarded as outlier.
% It is expected to be a p*1 array, where p is the number of channels.
% Each entry of outlierThr is expected to be [0, 0.01]
% The default value is 0.001 for each entry in OUTLIERTHR.
%
% The optional input position_coeff determines the ratio of the Euclidean
% between the center and each pixel.The default is 1.
%
% Output:
% The optional output myDataNor is a n*m array, that is ready for PCA
% analysis and clustering
% n: the biofilm pixel number; m: number of channels plus 1
% The end column is the distance towrads the biofilm center
% offer the position-information.
% Outliers are firstly removed then normalize to (0,1)
%
% The optional output myData is the rawdata (without normalization) of
% myDataNor
%
% If the output is not assigned to any variable, the function shows the
% histogram of all channels after normalization (myDataNor)
%
% offer the position-information.
%
% Write by Ping Shen and Modified by Tianmin Wang @ Tsinghua University
% Version 0.0.2. Created on Jul 27, 2019. Last modified on Aug 11, 2019.

% the optional input
argin = inputParser;
addParameter(argin,'outlierThr',0.001*ones(length(channelsOfInterest), 1))
addParameter(argin,'position_coeff',1)
parse(argin,varargin{:})
outlierThr = argin.Results.outlierThr;
position_coeff = argin.Results.position_coeff;

if length(outlierThr) ~= length(channelsOfInterest)
    error('Inconsistent number of entries for outlierThr!');
end

for m = 1:numel(channelsOfInterest)
    this_outlierThr = outlierThr(m);
    if this_outlierThr > 0.01
        error('incorrect outlier threshold setting!')
    elseif outlierThr < 0
        this_error('incorrect outlier threshold setting!')
    end
end

if ~iscell(M)
    error('M must be a cell.')
end
    
if length(M) ~= numel(channelsOfInterest)
    error('The number of channels in data is not consistent')
end

colonyPixels = length( find(colony) ); % number of colony pixels
myData = zeros(colonyPixels, (numel(channelsOfInterest) + 1));
    
% get the colony info for each channel
for m = 1:numel(channelsOfInterest)
    this_outlierThr = outlierThr(m); % threshold for this channel
    % data without normalization
    myData(:, m) = double( M{m}(colony));
    
    % remove the outlier pixel with extraordinary high FL
    Ms = sort(myData(:, m),'descend');   % Sort Descending
    outliers = Ms(1:ceil(length(Ms)*this_outlierThr));   % Top outliers
    upperLimit = outliers(end);
    myData(myData(:, m) > upperLimit, m) = upperLimit;
end

% add the radius info to myData
myData(:, end) = distance2Edge(colony);

% normalization
myDataNor = normalize(myData, 'range');
myDataNor(:, end) = myDataNor(:, end)*position_coeff; % adjust the wieght of position

% output
if nargout == 0
    % present the histogram of each channel
    colors = struct;
    colors.CFP = 'b';
    colors.GFP = 'g';
    colors.mCherry = 'm';
    colors.APC = 'r';
    edges = [0:0.02:1];
    clf;
    hold on;
    for m = 1:numel(channelsOfInterest)
        channel = channelsOfInterest{m};
        color = colors.(channel); 
        histogram(myDataNor(:, m), edges, ...
                  'DisplayStyle','stairs', 'EdgeColor',color);
    end
    hold off;
    alpha(0.3);
    xlabel('Fluorescence', 'Fontsize',18)
    ylabel('Number of cells', 'Fontsize',18)
    lgd = legend(channelsOfInterest, 'Location','northeast');
    lgd.FontSize = 14;
else
    varargout{1} = myDataNor;
    varargout{2} = myData;
end
