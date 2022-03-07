function labeledFACS = facsLabelAssign(dataImg, dataFACS, varargin)
% Assign each FACS cell with an index, according to its FL profile.
% Four steps in this fucntion:
% Extract N cells and N pixels from the rawdata in a random manner;
% Map the two extracted datasets according to optimal transport algorithm;
% In this way, label each FACS cell according to the label of its assigned
% image pixel.
% Optimize the categorizing by a 'the-most-distant-outlier-regrouping'
% strategy.
%
% Usage: labelFACS = facsLabelAssign(dataImg, dataFACS)
%        labelFACS = facsLabelAssign(..., 'maxSubset',1000)
%
% The input dataImg is a n*(p+1) matrix, rows specifying pixels and columns
% specify FL channels and label (last column).
%
% The input dataFACS is a m*p matrix, rows specifying cells and columns
% specify FL channels.
%
% Note that the order of FL channels should be the same for these two
% datasets.
%
% Optional input maxSubset specifies the max size of subset to split FACS
% dataset. The default is 1000. Bigger than 1000 will greatly slow down the
% calculation of optimal transport plan.
%
% Optional input OTsamplingFold specifies the size of the sampled imgData
% subset vs. size of the facsData when doing optimal transport.
% OTsamplingFold = 1 suggests strict abundance constraint for clustering.
% Bigger OTsamplingFold than 1 suggests relaxed abundance constrainst.
%
% The optional input compactThr determines the threshold defining whether a
% cluster is compact or not. This is used in the last optimization step.
% See outlierAssign function for details.
%
% The optional input maxIteration determines the maximum steps of
% iterations during the last optimization step.
%
% The optional input channels determines channels used to perform OT
% mapping. It is a 1*C array with min >= 1 and max <= size(dataImg, 2)-1.
% Default is 1:size(dataImg, 2) (All channels are used).
%
% Output:
% The output labeledFACS is a m*(p+1) array, adding the calculated label
% for each FACS cell as the last column.
%
% Dependency: OPmapping and outlierRegroup by Tianmin Wang
%
% Written by Tianmin Wang @ Tsinghua University
% Version 0.0.2. Created on Aug 02, 2019. Last modified on Jun 25, 2019.

argin = inputParser;
addOptional(argin,'maxSubset',1000)
addOptional(argin,'OTsamplingFold',3)
addOptional(argin,'compactThr',2)
addOptional(argin,'maxIteration',1000)
addOptional(argin,'channels',1:(size(dataImg,2)-1))
parse(argin,varargin{:})
maxSubset = argin.Results.maxSubset;
OTsamplingFold = argin.Results.OTsamplingFold;
compactThr = argin.Results.compactThr;
maxIteration = argin.Results.maxIteration;
channels = argin.Results.channels;

% check whether two datasets have the same amount of channels
[Npixels, ncImg] = size(dataImg);
[Ncells, ncFACS] = size(dataFACS);

if ncImg ~= ncFACS + 1
    error('inconsistent number of channels in two datasets!')
end

% check whether image data is large enough
if Npixels < maxSubset
    error('image data is not big enough for subsampling!')
end

% check whether OTsamplingFold is bigger or equal to 1
if OTsamplingFold < 1
    error('OTsamplingFold should be >= 1!')
end

% check whether channels is consistent with given data
if min(channels) < 1
    error('Channels must be a 1*C array specifying C channels used in OT mapping! (1=<C<size(dataImg,2))')
end

if max(channels) >= size(dataImg, 2)
    error('Channels must be a 1*C array specifying C channels used in OT mapping! (1=<C<size(dataImg,2))')
end

% shuffle FACS dataset
idxFACS = randperm(Ncells);
dataFACS_shuffle = dataFACS(idxFACS, :);

% final output
labeledFACS = [dataFACS_shuffle, zeros(Ncells, 1)];

% assign label to each FACS data subset
% split FACS data into subsets is due to the speed of optimal transport
% algorithm
for m = 1:ceil(Ncells/maxSubset)
    
    % split the FACS data
    if m*maxSubset <= Ncells
        FACSsubset = dataFACS_shuffle((m*maxSubset - maxSubset + 1):(m*maxSubset), :);
        subsetSize = maxSubset;
        lastCell = m*maxSubset;
    else
        FACSsubset = dataFACS_shuffle((m*maxSubset - maxSubset + 1):end, :);
        subsetSize = mod(Ncells, maxSubset);
        lastCell = Ncells;
    end
    startCell = m*maxSubset - maxSubset + 1;
    
    % randomly extract subsets from image data
    extracted_dataImg = datasample(dataImg, ceil(subsetSize*OTsamplingFold), ...
                                   'Replace',false);
    extracted_dataImg_ot = extracted_dataImg(:, 1:(end-1));
    extracted_imgLabel = extracted_dataImg(:, end);
    
    % optimal transport to assign labels to FACS subset
    [assign, cost] = OTmapping(extracted_dataImg_ot(:, channels), ...
                               FACSsubset(:, channels));
    
    % assign label to FACS cells
    labeledFACS(startCell:lastCell, end) = extracted_imgLabel(assign);
        
end

% reorder the data to make it consistent with input (delete shuffling)
for n = 1: ncFACS+1
    labeledFACS(idxFACS,n) = labeledFACS(1:end,n);
end