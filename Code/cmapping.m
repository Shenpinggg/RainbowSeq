function cmapping(cluster_index,colony,filename)
%
% This function is to map the clusters which is calculated by kmeans to the
% biofilm snapshot.
%
% Usage: cmapping(cluster_index,colony,filename)
%
% The input cluster_index should be a n*1 colomun vecter which is
% calculated by function kmeans.
% The input colony should be a logical matrix represent the biofilm pixels.
% The input filename is the integrate path including the format for storage
% such as 'E:\Data\D0500\mapping.png'.
% Written by Ping Shen @ Tsinghua University 
% Version 0.0.1. Created on Jul 15, 2019. Last modified on Aug 15, 2019.

% choose the colormap
figure
cmp = colormap(jet);
cmp = cmp(end:-1:1, :);
interval = floor(64/max(cluster_index));
color_map = cmp(interval*(1:max(cluster_index)),:);

% label the binary matrix colony with given cluster index
if  islogical(colony)
    idx_c = find(colony == 1);
    colony_new = zeros(size(colony));
    colony_new(idx_c) = cluster_index;
else
    error('The input colony should be a logical matrix')
end

% exhibit the picture of biofilm with cluster infromation and save the picture
new_image = label2rgb(colony_new,color_map);
imshow(new_image)
saveas(gcf,[filename '\mapping_to_biofilm.png'])
