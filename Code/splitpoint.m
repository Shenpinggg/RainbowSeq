function varargout = splitpoint(tree, myData_facs, myData_facs_nor, channelsOfInterest, varargin)

% This function is to replace the tree.CutPoint (which is the normalized data) by the raw FACS data.
% 
% Usage: [SplitPoint, boundTable] = splitpoint(tree, myData_facs, channelsOfInterest)
%        [SplitPoint, boundTable] = splitpoint(tree, myData_facs, channelsOfInterest, 'lbThr',0.01)
%
% Input:
% The input tree must be the decision tree calculated by function fitctree.
% 
% The input myData_facs is a n*m array, storing the rawdata of FACS.
% 
% The input myData_facs_nor is a n*m array, storing the data of FACS which has been normalized.
%
% The input channelsOfInterest is a cell array, specifying the channels of
% interest. Accepted channels include 'CFP', 'GFP', 'mCherry' and 'APC'
%
% The optional input lbThr is the lower boundary of each axis in the FACS
% space. The ubThr (upper boundary) is calculated by 1-lbThr. The default is
% 0.01.
%
% Output:
% The output SplitPoint is a n*1 vector represent the cut point to
% replace the tree.CutPoint.
%
% The input boundTable is a n*3 table storing the boundary values for each
% axis along the FACS space. The n is the number of clustering channels,
% thus all the same as the number of variables used to build the tree. See
% splitpoint.m for details.
%
% Written by Ping Shen @ Tsinghua University
% Version 0.0.3. Created on Aug 15, 2019. Last modified on Jun 25, 2020.

argin = inputParser;
addParameter(argin, 'lbThr', 0.01)
parse(argin, varargin{:});
lbThr = argin.Results.lbThr;
ubThr = 1 - lbThr;

if lbThr <= 0 || lbThr > 0.1
    error('incorrect lbThr setting!')
end

myData_facs_unique = struct;
myData_facs_nor_unique = struct;
myData_facs_nor_unique_sorted = struct;
myData_facs_unique_sorted = struct;

% prepare the facs dataset, perturbing the duplicate data points moderately
for n = 1:numel(channelsOfInterest)
    myData_facs_nor_unique.(channelsOfInterest{n}) = myData_facs_nor(:,n) + ...
                                                     (linspace(0, 0.00001, length(myData_facs_nor(:,n))))';
    myData_facs_nor_unique_sorted.(channelsOfInterest{n}) = sort(myData_facs_nor_unique.(channelsOfInterest{n}));                                            
    myData_facs_unique.(channelsOfInterest{n}) = myData_facs(:, n)+ (linspace(0, 0.001, length(myData_facs(:,n))))';
    myData_facs_unique_sorted.(channelsOfInterest{n}) = sort(myData_facs_unique.(channelsOfInterest{n}));
end 

% get the Node location in the tree's properties
Node_loc = find(tree.IsBranchNode == 1);
SplitPoint = zeros(numel(tree.IsBranchNode),1);

% use the interp1 to get the cutpoints in the FACS raw data
for n = 1:numel(Node_loc)
    x = Node_loc(n);
    SplitPoint(x,1) = interp1(myData_facs_nor_unique_sorted.(tree.CutPredictor{x}), myData_facs_unique_sorted.(tree.CutPredictor{x}), tree.CutPoint(x));
end

% store the upper and lower boundary of FACS space along each axis
boundTable = cell2table(channelsOfInterest', 'VariableNames',{'channel'});
boundTable.lb = zeros(length(channelsOfInterest), 1);
boundTable.ub = zeros(length(channelsOfInterest), 1);
[row, ~] = size(myData_facs);
for n = 1:numel(channelsOfInterest)
    channel = channelsOfInterest{n};
    lb_value = myData_facs_unique_sorted.(channel)(floor(row*lbThr));
    ub_value = myData_facs_unique_sorted.(channel)(ceil(row*ubThr));
    boundTable(n, 'lb') = {lb_value};
    boundTable(n, 'ub') = {ub_value};
end

% output
varargout{1} = SplitPoint;
varargout{2} = boundTable;
