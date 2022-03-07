function clusterInfo = calClusterInfo(facs_index, clustering_index, image_index, filename)

% This function is to calculate the purity, yield and abundance in the population of each cluster;
%
% Usage: clusterInfo = clusterInfo(facs_index, clustering_index, image_index, filename)
%
% Input
% The input facs_index should be a n*1 array with n entries, specifying the
% cluster number of each entry after OT mapping.
%
% The input clustering_index should be a n*1 array with n entries,
% specifying the cluster number of each entry after clustering using method like
% CART or MA-LS-chain.
%
% The input image_index should be a m*1 array with m entries, specifying
% the cluster number of each entry in the raw image data.
%
% Note that the cluster number can be positive integers or zero ( does not
% belong to any cluster).
%
% The input filename is the integrate path including the format for storage
% such as 'E:\Data\D0500'.\
%
% Output
% The output clusterInfo is a struct, with field 'ClusterX' specifying a
% table specifying yield and purity for each cluster.
%
% Written by Tianmin Wang @ Tsinghua University
% Version 0.0.1. Created on Oct 13, 2019. Last modified on Oct 13, 2019.

% check the consistency of row number of facs_index and clustering_index
[rowF, ~] = size(facs_index);
[rowC, ~] = size(clustering_index);
[rowI, ~] = size(image_index);

if rowF ~= rowC
    error('inconsistent number of entries in facs_index and clustering_index')
end

% check the consistency of cluster number
cNumF = max(facs_index);
cNumC = max(clustering_index);
cNumI = max(image_index);

if cNumF ~= cNumC || cNumF ~= cNumI || cNumC ~= cNumI
    error('inconsistent number of clusters!')
else
    cNum = cNumF;
end

% plot the relative abundance of each cluster
ratio_img = zeros(cNum, 1);
ratio_facs = zeros(cNum, 1);
ratio_clustering = zeros(cNum, 1);
xticklabel_name = cell(cNum, 1);
for m = 1:cNum
    ratio_img(m) = numel(find(image_index == m))/rowI;
    ratio_facs(m) = numel(find(facs_index == m))/rowF;
    ratio_clustering(m) = numel(find(clustering_index == m))/rowC;
    xticklabel_name{m} = ['Cluster ' num2str(m)];
end
ratio_img_facs = [ratio_img ratio_facs ratio_clustering];
figure
bar(ratio_img_facs)
legend('IMAGE', 'FACS', 'Final')
title('The ratio of each cluster', 'FontSize', 20)
ylabel('Probability (%)','Fontsize',18)
set(gca,'Xticklabel',xticklabel_name,'FontSize',12)
xtickangle(90)
saveas(gcf, [filename '\abundance.png'])

% calculate the purity and yield of each cluster
clusterInfo = struct;
for m = 1:cNum
    cluster = ['Cluster', num2str(m)];
    clusterInfo.(cluster) = array2table(zeros(1,2), 'VariableNames',{'purity', 'yield'});
    purity = length(find(facs_index == m & clustering_index == m)) / length(find(clustering_index == m));
    yield = length(find(facs_index == m & clustering_index == m)) / length(find(facs_index == m));
    clusterInfo.(cluster)(1, 'purity') = {purity};
    clusterInfo.(cluster)(1, 'yield') = {yield};
end