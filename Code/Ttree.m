function tree = Ttree(mydata, index, channelsOfInterest, NumSplits, varargin)
%
% This function is to build up a binary tree according to the given data 
% and index.
% Usage:
%       tree = Ttree(mydata, index, channelsOfInterest, NumSplits)
%       tree = Ttree(mydata, index, channelsOfInterest, NumSplits, 'NumKfold', ...)
% 
% Input:
% The input mydata is the normalized data contained fluorescense intensity
% from each channels.
% The input index contains the result of clustering.
% The input channelsOfInterest is a cell array, specifying the channels of
% interest. Accepted channels include 'CFP', 'GFP', 'mCherry' and 'APC'.
% The input NumSplits should be a scalar represent the tree's number for
% spliting, recommended number of clusters plus 20.
% The optional input NumKfold determine the group numbers for K-fold 
% Cross-validate, default is 10.
% 
% Output:
% The noly output tree is the binary tree calculated by function fitctree.
%
% Written by Ping Shen @ Tsinghua University
% Version 0.0.2. Created on Aug 19, 2019. Last modified on Jun 25, 2020.
 

argin = inputParser;
addParameter(argin, 'NumKfold', 10);
parse(argin, varargin{:});
NumKfold = argin.Results.NumKfold;

% load data for cartTree
data_table = array2table(mydata, 'VariableNames', channelsOfInterest);
group = index;
length_data = size(data_table, 1);

% divide data into 2 parts, one for CartTree training, another for cross
% validation
idx = crossvalind('Kfold', length_data, NumKfold);
i = 1;
test = (idx == i);
train = ~test;
data_table_train = data_table(train, :);
group_train = group(train, :);
data_table_test = data_table(test, :);
group_test = group(test, :);

% Training CartTree
tree = fitctree(data_table_train, group_train, 'MaxNumSplits', NumSplits);
view(tree, 'Mode', 'graph')

% Predict the precision of CART tree based on data for validation
label = predict(tree, data_table_test);
a = find(label == group_test);
precision = length(a) / length(group_test);
sprintf(['The precision of the CART tree is ' , num2str(precision)])