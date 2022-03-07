function pca_plot(score, coeff, channelsOfInterest, filename, cluster_index)
%
% pca_plot(score, coeff, channelsOfInterest, filename, cluster_index)
% 
% This function can plot the pca associated data.
% The input score is principal component score which is calculated by
% function pca.If the channelsOfInterest contain only 2 elements(channls), 
% the input score can be replaced by input mydataNor which is a n*p matrix 
% calculated by imgMat. 
%
% The input coeff represent the principal component coefficients calculated
% by function pca from the same data.If the channelsOfInterest only contain
% 2 elements, the input coeff would not be used in the further analysis and
% user should set a random object,such as [].
%
% The input filename is the integrate path including the format for storage
% such as 'E:\Data\D0500\mapping.fig'.
%
% The input cluster_index ccontain cluster information calculated from
% kmeans.
%
% Written by Ping Shen @ Tsinghua University 
% Version 0.0.1. Created on Jul 15, 2019. Last modified on Aug 15, 2019.

[row,col] = size(coeff);
new_axis = [zeros(1, col); coeff];
nchannels = numel(channelsOfInterest);
% assign color
colors = struct;
colors.CFP = 'b';
colors.GFP = 'g';
colors.mCherry = 'm';
colors.APC = 'r';

figure
if numel(channelsOfInterest) == 2
   cmp = colormap(jet);
   cmp = cmp(end:-1:1, :);
   interval = floor(64/max(cluster_index));
   hold on
   for n = 1:max(cluster_index)
       color_map(n,:) = cmp(interval*n,:);
       histogram2(score(cluster_index == n,1), score(cluster_index == n,2),'DisplayStyle','bar3','FaceColor',color_map(n,:),'EdgeAlpha',0);
       cluster_name{n} = sprintf('cluster %d',n);
   end
   xlabel(channelsOfInterest{1})
   ylabel(channelsOfInterest{2})
   zlabel('Number of pixels')
   title('Clusters');
   legend(cluster_name);
   alpha(0.5);
   grid on
   savefig([filename '\clusters.fig'])
else % pca fluorescence plot with cluster information
    clf;
    histogram2(score(:,1), score(:,2), 'FaceColor','flat'); xlabel('PCA 1'); ylabel('PCA 2'); zlabel('Number of pixels');
    colorbar
    title('PCA (fluorescence)');
    saveas(gcf,[filename '\pca_fluo.png'])
    
    clf;
    hold on
    histogram2(score(:,1), score(:,2),'DisplayStyle','tile', 'FaceColor','flat'); xlabel('PCA 1'); ylabel('PCA 2'); zlabel('Number of pixels');
    for i = 1:nchannels
        color = colors.(channelsOfInterest{i});
        plot(new_axis([1 i+1],1), new_axis([1 i+1],2), 'LineWidth',5, 'Color',color)
    end
    colorbar
    grid on
    hold off
    saveas(gcf, [filename '\pcatile_fluo.png'])
    
    % pca cluster plot
    clf;
    cmp = colormap(jet);
    cmp = cmp(end:-1:1, :);
    interval = floor(64/max(cluster_index));
    hold on
    for n = 1:max(cluster_index)
        color_map(n,:) = cmp(interval*n,:);
        histogram2(score(cluster_index == n,1), score(cluster_index == n,2),'DisplayStyle','bar3','FaceColor',color_map(n,:),'EdgeAlpha',0);
        cluster_name{n} = sprintf('cluster %d',n);
    end
    alpha(0.5);
    for i = 1:nchannels
        color = colors.(channelsOfInterest{i});
        plot(new_axis([1 i+1],1), new_axis([1 i+1],2), 'LineWidth',5, 'Color',color)
    end
    hold off
    title('PCA (cluster)');
    legend(cluster_name);
    xlabel('PCA 1'); ylabel('PCA 2'); zlabel('Number of pixels');
    grid on
    savefig([filename '\pca_cluster.fig'])
end
clf;
