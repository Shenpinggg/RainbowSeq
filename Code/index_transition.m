function index_sorted = index_transition(cluster_index, reference, varargin)

% Usage: index_sorted = index_transition(cluster_index, reference)
%        index_sorted = index_transition(cluster_index, reference, ...
%                                        'Direction', ...)
% The index_transition is to resort the cluster_index calculated by kmeans
% based on the given reference.
% 
% The input cluster_index should be a column vector contained the cluster
% infromation.
% The input reference should also be a column vector as the basal to adjust
% the order of cluster index.
% 
% The optional input Direction determines the way of sorting.(Default is 0)
%              1 ---- 'ascent'
%              0 ---- 'descent'
%
% The output index_sorted is the new cluster index that has been resorted.
%
% Written by Ping Shen @ Tsinghua University
% Version 0.0.1. Created on Aug 15, 2019.


argin = inputParser;
addParamValue(argin,'Direction',0)
parse(argin,varargin{:})
Direction = argin.Results.Direction;


max_value = max(cluster_index);
order = [1:max_value];

% Sort the intensity
for n = 1:max_value
    value(n,1) = mean(reference(cluster_index == n));
end

switch Direction
    case 1
        [~,I_order] = sort(value);
    case 0
        [~,I_order] = sort(value,'descend');
end
% Store the original cluster index
for n = 1:max_value
    list_name{n} = sprintf('idx%d',n);
    Cluster.(list_name{n}) = find(cluster_index == n);
end

% Change the order 
for n = 1:max_value
    cluster_index(Cluster.(list_name{I_order(n)})) = order(n);
end

index_sorted = cluster_index;


