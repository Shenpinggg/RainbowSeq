function varargout = facsMat(biofilmFACSdata, ncFACSdata, channelsOfInterest, varargin)

% This function is to process biofilm facs data to a n*m array, ready for
% PCA analysis and clustering;
% n: bacterial cells; m: number of channels
%
% Usage: [myDataNor, myData] = facsMat(biofilmFACSdata, ncFACSdata, channelsOfInterest)
%        facsMat(biofilmFACSdata, ncFACSdata, channelsOfInterest)
%        [myDataNor, myData] = facsMat(..., 'ncThr',[1,1,1])
%        [myDataNor, myData] = facsMat(..., 'dataThr',9000)
%        [myDataNor, myData] = facsMat(..., 'outlierThr',[0.001,0.01,0.01])
%        [myDataNor, myData] = facsMat(..., 'Iflog',1)
%        [myDataNor, myData] = facsMat(..., 'FSCnor',false)
%
% Input
% The input FACSdata should be a standard .fcs file with specified path.
% One such data corresponds to the biofilm (biofilmFACSdata); 
% another is the negative control(without dye staining) (ncFACSdata)
%
% The input channelsOfInterest is a cell array, specifying the channels of
% interest. Accepted channels include 'CFP', 'GFP', 'mCherry' and 'APC'
%
% The optional input ncThr determines the fold of nc median used to
% normalize the data. For a particular channel, say, GFP, we have median of
% nc cells, say, medianNC. The biofilm data is processed by the equation:
% biofilm(biofilm('GFP') < medianNC*ncThr) = medianNC*ncThr
% It is expected to be a p*1 array, where p is the number of channels.
% Each entry in ncThr is expected to be > 1
% The default value is 1 for each entry of NCTHR.
%
% The optional input dataThr determines largest number of data entries to
% be processed. If raw data has more entries than dataThr, random
% subsampling will be performed. Default: 9000.
%
% The optional input outlierThr determines largest proportion of data in
% one channel, that is regarded as outlier.
% It is expected to be a p*1 array, where p is the number of channels.
% Each entry of outlierThr is expected to be [0, 0.01]
% The default value is 0.001 for each entry of OUTLIERTHR.
% 
% The optional input Iflog determines whether pre-treat the raw facs data
% with function log10(), the 0 means non pretreatment, the 1 means
% log10() pretreatment.(The default is 1)
% 
% The optional input FSCnor specifies whether fluorescence data need to be
% normalized by FSC. The equation is fluorescence/FSC*median(FSC). Default
% is false.
%
% Output
% The optional output myDataNor is a n*m array, that is log10 processed and
% normalized, ready for PCA analysis and clustering;
% n: bacterial cells; m: number of channels
% the data is firstly projected onto log10 space and then normalize to (0,1)
%
% The optional output myData is a n*m array, storing the rawdata of FACS;
% n: bacterial cells; m: number of channels
%
% If the output is not assigned to any variable, the function shows the
% histogram of all channels after normalization (myDataNor)
%
% Dependency: facsRead by Tianmin Wang & fca_readfcs by Laszlo Balkay
%
% Written by Tianmin Wang @ Tsinghua University
% Version 0.0.5. Created on Jul 27, 2019. Last modified on Jun 25, 2020.

% optional
argin = inputParser;
addParameter(argin,'ncThr',ones(length(channelsOfInterest)))
addParameter(argin,'outlierThr',0.001*ones(length(channelsOfInterest), 1))
addParameter(argin,'dataThr',9000)
addParameter(argin,'Iflog',1)
addParameter(argin,'FSCnor',false)
parse(argin,varargin{:})
ncThr = argin.Results.ncThr;
outlierThr = argin.Results.outlierThr;
dataThr = argin.Results.dataThr;
Iflog = argin.Results.Iflog;
FSCnor = argin.Results.FSCnor;

if length(ncThr) ~= length(channelsOfInterest)
    error('Inconsistent number of entries for ncThr!');
end

if length(outlierThr) ~= length(channelsOfInterest)
    error('Inconsistent number of entries for outlierThr!');
end

for m = 1:numel(channelsOfInterest)
    this_ncThr = ncThr(m);
    this_outlierThr = outlierThr(m);
    if this_ncThr < 1
        error('incorrect control threshold setting!')
    end
    if this_outlierThr > 0.01
        error('incorrect outlier threshold setting!')
    elseif outlierThr < 0
        error('incorrect outlier threshold setting!')
    end
end

% read the raw data
biofilmData = facsRead(biofilmFACSdata);
ncData = facsRead(ncFACSdata);

% check the channel consistency
mychannels = biofilmData.Properties.VariableNames;
for m = 1:numel(channelsOfInterest)
    channel = channelsOfInterest{m};
    if ~contains(mychannels, channel)
        error(('incorrect channel name: %s'), channel);
    end
end

% get a cell array containing all channels of interest with FSC as the last
% item
channelsPlusFSC = channelsOfInterest;
channelsPlusFSC{end + 1} = 'FSC';

myBiofilm = biofilmData(:, channelsPlusFSC);
myNC = ncData(:, channelsPlusFSC);

myData = table2array(myBiofilm);

% get subset of the data
dataThr = min(dataThr, length(myData));
myData = datasample(myData, dataThr, 'Replace',false);

% normalization of myData by FSC-A
myDataNor = myData;
if FSCnor
    FSCcolumn = myData(:, end)./ median( myData(:, end) ); % normalize vector
    myDataNor = myData./FSCcolumn;
end
% delete the FSC column
myDataNor = myDataNor(:, 1:end-1);
myData = myData(:, 1:end-1);

% normalization of myNC by FSC-A
if FSCnor
    FSCcolumn = table2array(myNC(:, 'FSC'));
    FSCcolumn = FSCcolumn ./ median(FSCcolumn); % normalize vector
    for c = 1:length(channelsOfInterest)
        channel = channelsOfInterest{c};
        myNC(:, channel) = array2table( table2array(myNC(:, channel)) ...
                                        ./ FSCcolumn);
    end
end
% delete the FSC column
myNC = myNC(:, channelsOfInterest);

% normalize the data by minus the negative control
for m = 1:numel(channelsOfInterest)
    this_ncThr = ncThr(m);   % threshold for this channel
    channel = channelsOfInterest{m};
    ncMedian = max(median(table2array(myNC(:, channel))), 1/this_ncThr); % make the minimum as 1, thus 0 after log10 processing
    ncMedian = ncMedian*this_ncThr;
    myDataNor(myDataNor(:, m) < ncMedian, m) = ncMedian;
end

% replace outliers in the data of each channel
for m = 1:numel(channelsOfInterest)
    this_outlierThr = outlierThr(m);  % threshold for this channel
    Ms = sort(myDataNor(:, m),'descend');                % Sort Descending
    outliers = Ms(1:ceil(length(Ms)*this_outlierThr));    % Top outliers
    upperLimit = outliers(end);
    myDataNor(myDataNor(:, m) > upperLimit, m) = upperLimit;
end

% log10 and normalization
if Iflog
    myDataNor = log10(myDataNor);
end
    
% project to [0, 1]
myDataNor = normalize(myDataNor, 'range');

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
        % histogram(myDataNor(:, m), 'FaceColor',color);
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


