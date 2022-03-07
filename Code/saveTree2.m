function varargout = saveTree2(FACS_idx, myData_facs, tree, SplitPoint, boundTable, channelsOfInterest, filename, varargin)

% This function is to conclude the tree's properties and save the split
% pathway of each leaf node belonging to one particular cluster. Note that
% one cluster has >= 1 leaf node(s).
% 
% Usage: saveTree2(FACS_idx, myData_facs, tree, SplitPoint, boundTable, channelsOfInterest, filename)
%        [final_Decision, raw_tree, clusterBoundInfo, xmlCluster] = saveTree2(...)
%        [final_Decision, raw_tree, clusterBoundInfo, xmlCluster] = saveTree2(..., 'ViewTree', 0)
%
% Input
% 
% The input FACS_idx specifies the group of each cell assigned by OT
% mapping
% 
% The input myData_facs is a n*m array, storing the rawdata of FACS. The
% order of columns corresponds to the input variable "channelsOfInterest"
% 
% The input tree must be the decision tree calculated by function fitctree.
% 
% The input SplitPoint should be a n*1 vector represented the cut points to
% replace the tree.CutPoint.
% 
% The input boundTable is a n*3 table storing the boundary values for each
% axis along the FACS space. The n is the number of clustering channels,
% thus all the same as the number of variables used to build the tree. See
% splitpoint.m for details.
%
% The input channelsOfInterest is a cell array, specifying the channels of
% interest. Accepted channels include 'CFP', 'GFP', 'mCherry' and 'APC'
%
% The input filename is the integrate path including the format for storage
% such as 'E:\Data\D0500'.\
%
% The optional input ViewTree is to choose whether exhibit the view-graph
% and view-text of the given tree.
%
% Output
% The clusterBoundInfo is a struct with 'cluster X' as field, pointing to a cell array containing >= 1 table(s).
% Number of items in the cell array corresponds to the number of leaf nodes
% calculated by CART tree belonging to this very cluster. In the cell array,
% each table has n*3 elements, where n is number of channels used in
% clustering. The first column is the channel, the second is the lower
% boundary and the thrid is the upper boundary.
%
% The xmlCluster is a struct mimicking the xml file, describing two 4-point
% polygon gates used to obtain each leaf node belonging to one cluster. In particular, it is a
% hierarchical struct with structures like the following.
% ClusterX -> a cell array, in this cell array, we have a hierarchical
% struct, with two fields: gates and properties.
% Firstly, gates{gate1, gate2, ...} (gate is a struct) -> xparam, yparam,
% used (whether it will be used to write in the final .xml file), 
% points{point1, point2, point3, point4} (point is a struct) -> x, y
% Thus, you can access the first point of gate1 of cluster 1 using:
% xmlCluster.Cluster1.gates{1}.points{1}.x(or y)
% Note that here is just the gates calculated from CART tree. It can be
% further refined by additional processing.
% Secondly, properties is a struct, pointing to three fields, purity,
% relativeYield and absoluteYield. Purity simply represents the purity for
% this leaf node; relativeYield represents the yield of this leaf node
% relative to the cells in the relevant cluster; absoluteYield represents
% the yield of this leaf node relative to all the cells.
% 
% Written by Ping Shen and Tianmin Wang @ Tsinghua University
% Version 0.3.1. Created on Aug 15, 2019. Last modified on Jun 30, 2020.

argin = inputParser;
addParameter(argin, 'ViewTree', 0)
parse(argin, varargin{:});
ViewTree = argin.Results.ViewTree;

% key propersities of the tree for getting the integrated split-path
raw_tree = table(tree.IsBranchNode,tree.CutPoint,tree.CutPredictor,tree.NodeClass,tree.Children,tree.Parent,...
    'VariableNames',{'IsBranchNode','CutPoint','CutPredictor','NodeClass','Children','Parent'});

switch ViewTree
    case 1
    view(tree)
    view(tree, 'Mode', 'Graph')
    case 0
end

Leaf_loc = find(tree.IsBranchNode == 0);
Node_loc = find(tree.IsBranchNode == 1);
Decision_tree = cell(numel(tree.IsBranchNode),2);

% create a cell contained different cut decisions
for n = 1:numel(tree.IsBranchNode)
    if tree.IsBranchNode(n)
       Decision_tree{n,1} = [tree.CutPredictor{n}, ' < ', num2str(SplitPoint(n))];
       Decision_tree{n,2} = [tree.CutPredictor{n}, ' >= ', num2str(SplitPoint(n))];
    else
       Decision_tree{n,1} = '';
       Decision_tree{n,2} = '';
    end
end

LeafClass = cell(numel(Leaf_loc),1);
SplitPath = cell(numel(Leaf_loc),numel(Node_loc));

% combine the pathes to a integrate split-path
for n = 1:numel(Leaf_loc)
    x = Leaf_loc(n);
    LeafClass{n} = ['Cluster', tree.NodeClass{x}]; 
    m = 1;
    while find(tree.Children == x)
    [x, y] = find(tree.Children == x);
    SplitPath{n,m} = Decision_tree{x,y};
    m = m + 1;
    end
end

final_Decision = table(LeafClass, SplitPath);
writetable(final_Decision, [filename '\treeDecision.csv']);

% prepare boundary information for each cluster obtained from CART tree,
% used in further refinement using MS_LS-chain.

% matrix of boundary along each axis of all cells
lb_array = zeros(length(channelsOfInterest), 1);
ub_array = zeros(length(channelsOfInterest), 1);
for n = 1:numel(channelsOfInterest)
    channel = channelsOfInterest{n};
    lb_array(n,1) = table2array(boundTable(strcmp(boundTable.channel, channel), 'lb'));
    ub_array(n,1) = table2array(boundTable(strcmp(boundTable.channel, channel), 'ub'));
end

% record the number of clusters and the initialization of clusterBoundInfo,
% clusterChangedChannel
clusterBoundInfo = struct;
all_clusters = {};
updated_LN_of_each_cluster = struct; % number of LN (leaf node) updated for each cluster
clusterActiveChannel = struct; % record about the active channel(s) for each leaf node, storing all unique active channels
for n = 1:numel(Leaf_loc)
    x = Leaf_loc(n);
    cluster = ['Cluster', tree.NodeClass{x}];
    % record the number of clusters and the initialization of clusterBoundInfo
    if isempty(find(strcmp(cluster, all_clusters), 1))
        all_clusters{end + 1} = cluster;
        clusterBoundInfo.(cluster) = {};
        clusterActiveChannel.(cluster) = {};
        updated_LN_of_each_cluster.(cluster) = 0;
    end
    % initialization of each leaf node belonging to one corresponding cluster in clusterBoundInfo
    clusterBoundInfo.(cluster){end + 1} = cell2table(channelsOfInterest', 'VariableNames',{'channel'});
    clusterBoundInfo.(cluster){end}.lb = lb_array;
    clusterBoundInfo.(cluster){end}.ub = ub_array;
    clusterActiveChannel.(cluster){end + 1} = {}; % initialize clusterActiveChannel
end

% replace the raw boundary (all cells) by CART boundary for each leaf node belonging to one particular cluster
for n = 1:numel(Leaf_loc)
    x = Leaf_loc(n);
    cluster = ['Cluster', tree.NodeClass{x}];
    updated_LN_of_each_cluster.(cluster) = updated_LN_of_each_cluster.(cluster) + 1; % the leaf node that will be updated
    LN_to_be_updated = updated_LN_of_each_cluster.(cluster); % this leaf node will be updated in clusterBoundInfo
    allBranch = final_Decision(strcmp(final_Decision.LeafClass, cluster), 'SplitPath'); % inner nodes for all leaf nodes belonging to this cluster
    % search the split path to update the boundary
    allBranch = table2cell(allBranch);
    thisBranch = allBranch{LN_to_be_updated}; % inner nodes for this leaf node to be processed
    for m = 1:numel(thisBranch)
        oneBranch = thisBranch{m};
        if ~isempty(oneBranch)
            if contains(oneBranch, ' < ')
                channel = extractBefore(oneBranch, ' < ');
                channel_ub = str2double(extractAfter(oneBranch, ' < '));
                channel_to_be_updated = strcmp(clusterBoundInfo.(cluster){LN_to_be_updated}.channel, channel);
                if channel_ub < table2array(clusterBoundInfo.(cluster){LN_to_be_updated}(channel_to_be_updated, 'ub')) % update only when the upper bound of this inner node is lower than the record
                    clusterBoundInfo.(cluster){LN_to_be_updated}(channel_to_be_updated, 'ub') = {channel_ub};
                end
                if isempty(find(strcmp(channel, clusterActiveChannel.(cluster){LN_to_be_updated}), 1)) % store the active channel for this leaf node
                    clusterActiveChannel.(cluster){LN_to_be_updated}{end + 1} = channel;
                end
            elseif contains(oneBranch, ' >= ')
                channel = extractBefore(oneBranch, ' >= ');
                channel_lb = str2double(extractAfter(oneBranch, ' >= '));
                channel_to_be_updated = strcmp(clusterBoundInfo.(cluster){LN_to_be_updated}.channel, channel);
                if channel_lb > table2array(clusterBoundInfo.(cluster){LN_to_be_updated}(channel_to_be_updated, 'lb')) % update only when the lower bound of this inner node is higher than the record
                    clusterBoundInfo.(cluster){LN_to_be_updated}(channel_to_be_updated, 'lb') = {channel_lb};
                end
                if isempty(find(strcmp(channel, clusterActiveChannel.(cluster){LN_to_be_updated}), 1)) % store the active channel for this leaf node
                    clusterActiveChannel.(cluster){LN_to_be_updated}{end + 1} = channel;
                end
            end
        end
    end
end

% sort active channels for each leaf node according to the
% channelsOfInterest variable
for k = 1:numel(all_clusters)
    cluster = all_clusters{k};
    for h = 1:numel(clusterActiveChannel.(cluster)) % update each leaf node belonging to this cluster
        [~, idx] = ismember(clusterActiveChannel.(cluster){h}, channelsOfInterest);
        ordered_idx = sort(idx);
        clusterActiveChannel.(cluster){h} = channelsOfInterest(ordered_idx);
    end
end

% the initialization of xmlCluster and re-initialization of updated_LN_of_each_cluster
xmlCluster = struct;
updated_LN_of_each_cluster = struct; % number of LN (leaf node) updated for each cluster
all_clusters = {};
for n = 1:numel(Leaf_loc)
    x = Leaf_loc(n);
    cluster = ['Cluster', tree.NodeClass{x}];
    % record the number of clusters and the initialization of clusterBoundInfo
    if isempty(find(strcmp(cluster, all_clusters), 1))
        all_clusters{end + 1} = cluster;
        xmlCluster.(cluster) = {};
        updated_LN_of_each_cluster.(cluster) = 0;
    end
    % initialization of each leaf node belonging to one corresponding cluster in clusterBoundInfo
    xmlCluster.(cluster){end + 1} = struct;
end

% update each leaf node in xmlCluster
xmlCluster = struct;
for n = 1:numel(Leaf_loc)
    x = Leaf_loc(n);
    cluster = ['Cluster', tree.NodeClass{x}];
    updated_LN_of_each_cluster.(cluster) = updated_LN_of_each_cluster.(cluster) + 1; % the leaf node that will be updated
    LN_to_be_updated = updated_LN_of_each_cluster.(cluster); % this leaf node will be updated in clusterBoundInfo
    
    % update the gates fields of each leaf node
    xmlCluster.(cluster){LN_to_be_updated}.gates = {};
    used_channels = {};
    for m = 1:numel(channelsOfInterest)
        channelx = channelsOfInterest{m};
        for k = 1:numel(channelsOfInterest)
            channely = channelsOfInterest{k};
            if ~strcmp(channelx, channely) && isempty(find(contains(used_channels, channely), 1)) % combine two channels to have a gate, Cnk
                xmlCluster.(cluster){LN_to_be_updated}.gates{end + 1} = struct; % new gate for this lead node
                channelx_to_be_updated = strcmp(clusterBoundInfo.(cluster){LN_to_be_updated}.channel, channelx); % row number of channelx in the table
                channely_to_be_updated = strcmp(clusterBoundInfo.(cluster){LN_to_be_updated}.channel, channely); % row number of channely in the table
                xmlCluster.(cluster){LN_to_be_updated}.gates{end}.xparam = channelx;
                xmlCluster.(cluster){LN_to_be_updated}.gates{end}.yparam = channely;
                xmlCluster.(cluster){LN_to_be_updated}.gates{end}.active = false;
                xmlCluster.(cluster){LN_to_be_updated}.gates{end}.points = {struct, struct, struct, struct};
                xmlCluster.(cluster){LN_to_be_updated}.gates{end}.points{1}.x = table2array(clusterBoundInfo.(cluster){LN_to_be_updated}(channelx_to_be_updated, 'lb'));
                xmlCluster.(cluster){LN_to_be_updated}.gates{end}.points{1}.y = table2array(clusterBoundInfo.(cluster){LN_to_be_updated}(channely_to_be_updated, 'ub'));
                xmlCluster.(cluster){LN_to_be_updated}.gates{end}.points{2}.x = table2array(clusterBoundInfo.(cluster){LN_to_be_updated}(channelx_to_be_updated, 'ub'));
                xmlCluster.(cluster){LN_to_be_updated}.gates{end}.points{2}.y = table2array(clusterBoundInfo.(cluster){LN_to_be_updated}(channely_to_be_updated, 'ub'));
                xmlCluster.(cluster){LN_to_be_updated}.gates{end}.points{3}.x = table2array(clusterBoundInfo.(cluster){LN_to_be_updated}(channelx_to_be_updated, 'ub'));
                xmlCluster.(cluster){LN_to_be_updated}.gates{end}.points{3}.y = table2array(clusterBoundInfo.(cluster){LN_to_be_updated}(channely_to_be_updated, 'lb'));
                xmlCluster.(cluster){LN_to_be_updated}.gates{end}.points{4}.x = table2array(clusterBoundInfo.(cluster){LN_to_be_updated}(channelx_to_be_updated, 'lb'));
                xmlCluster.(cluster){LN_to_be_updated}.gates{end}.points{4}.y = table2array(clusterBoundInfo.(cluster){LN_to_be_updated}(channely_to_be_updated, 'lb'));
            end
        end
        used_channels{end + 1} = channelx;
    end
    
    % update the properties field of each leaf node
    xmlCluster.(cluster){LN_to_be_updated}.properties = struct;
    [cellNum, ~] = size(myData_facs);
    record = true(cellNum, 1); % record about whether each cell belongs to this leaf node or not
    % find the cells belonging to this leaf node
    for m = 1:numel(channelsOfInterest)
        channel = channelsOfInterest{m};
        channelNum = strcmp(clusterBoundInfo.(cluster){LN_to_be_updated}.channel, channel); % the number of channel in this table
        channel_lb = table2array(clusterBoundInfo.(cluster){LN_to_be_updated}(channelNum, 'lb'));
        channel_ub = table2array(clusterBoundInfo.(cluster){LN_to_be_updated}(channelNum, 'ub'));
        record = record & ( (myData_facs(:, m) >= channel_lb) & (myData_facs(:, m) <= channel_ub));
    end
    % calculate purity and yield
    clusterOfThisLN = str2num(tree.NodeClass{x});
    purity = length(find(FACS_idx == clusterOfThisLN & record)) / length(find(record));
    relativeYield = length(find(FACS_idx == clusterOfThisLN & record)) / length(find(FACS_idx == clusterOfThisLN));
    absoluteYield = length(find(FACS_idx == clusterOfThisLN & record)) / cellNum;
    xmlCluster.(cluster){LN_to_be_updated}.properties.purity = purity;
    xmlCluster.(cluster){LN_to_be_updated}.properties.relativeYield = relativeYield;
    xmlCluster.(cluster){LN_to_be_updated}.properties.absoluteYield = absoluteYield;
    
end

% update the active gate(s) in xmlCluster
% the basic strategy is to find the minimal number of polygon gates (in 2D
% fluorescent space), whose channels cover all active channels for the
% relevant leaf node stored in clusterActiveChannel;
for k = 1:numel(all_clusters)
    cluster = all_clusters{k};
    for h = 1:numel(xmlCluster.(cluster)) % update each leaf node belonging to this cluster
        numActiveChannels = numel(clusterActiveChannel.(cluster){h});
        if numActiveChannels == 0
            error(['No active channel for ', cluster, ' leaf node ', num2str(h), '!'])
        elseif numActiveChannels == 1 % only one active channel, find another different channel to form the polygon gate
            if ~strcmp(channelsOfInterest{1}, clusterActiveChannel.(cluster){h}{1})
                channelx = channelsOfInterest{1};
                channely = clusterActiveChannel.(cluster){h}{1};
            else
                channelx = channelsOfInterest{1};
                channely = channelsOfInterest{2};
            end
            for n = 1:numel(xmlCluster.(cluster){h}.gates) % find the relevant gate for this leaf node with corresponding channels
                if strcmp(xmlCluster.(cluster){h}.gates{n}.xparam, channelx) && strcmp(xmlCluster.(cluster){h}.gates{n}.yparam, channely)
                    xmlCluster.(cluster){h}.gates{n}.active = true;
                    break
                end
            end
        else % at least two active channels, use polygon with active channel 1 + 2, 3 + 4, ... etc; the last single channel will form a polygon with the first channel in channelsOfInterest
            for m = 0:2:numActiveChannels
                if (numActiveChannels - m) >= 2
                    channelx = clusterActiveChannel.(cluster){h}{m + 1};
                    channely = clusterActiveChannel.(cluster){h}{m + 2};
                    for n = 1:numel(xmlCluster.(cluster){h}.gates) % find the relevant gate for this leaf node with corresponding channels
                        if strcmp(xmlCluster.(cluster){h}.gates{n}.xparam, channelx) && strcmp(xmlCluster.(cluster){h}.gates{n}.yparam, channely)
                            xmlCluster.(cluster){h}.gates{n}.active = true;
                            break
                        end
                    end
                elseif (numActiveChannels - m) == 1
                    channelx = channelsOfInterest{1};
                    channely = clusterActiveChannel.(cluster){h}{numActiveChannels};
                    for n = 1:numel(xmlCluster.(cluster){h}.gates) % find the relevant gate for this leaf node with corresponding channels
                        if strcmp(xmlCluster.(cluster){h}.gates{n}.xparam, channelx) && strcmp(xmlCluster.(cluster){h}.gates{n}.yparam, channely)
                            xmlCluster.(cluster){h}.gates{n}.active = true;
                            break
                        end
                    end
                elseif (numActiveChannels - m) == 0
                    break
                end
            end
        end
    end
end
              
% output
varargout{1} = final_Decision;
varargout{2} = raw_tree;
varargout{3} = clusterBoundInfo;
varargout{4} = xmlCluster;
