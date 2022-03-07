function cplot3_2(mydata, Orderofchannels, cluster_index, filename, varargin)
% 
% Display clustered data in 3-D space including unclustered sample.
% Usage: cplot3_2(mydata, Orderofchannels, cluster_index)
%
% The input mydata should be a n*m matrix, n represent the sample number
% and m represent the channel number. The channel order is suggested to
% CFP,GFP,mCherry...(the order is wavelength from short to long).
%
% The input Orderofchannels is a cell array contained 3 elements represent 
% each channel.The channel order is suggested to CFP,GFP,mCherry...
% (the order is wavelength from short to long).
%
% The input filename is the integrate path including the format for storage
% such as 'E:\Data\D0500\mapping.png'.
%
% The input cluster_index should be a n*1 colomun vecter which is
% predicted by CARTtree.The max index represents unclustered cells
% from prediction of CARTtree
% 
% The optional input show_unclustered means whether figure show the
% unclustered cells.The dafault value is 1 represent yes.
%
% Written by Ping Shen @ Tsinghua University 
% Version 0.0.1. Created on Jul 09, 2020.

argin = inputParser;
addOptional(argin,'show_unclustered',1)
parse(argin,varargin{:})
show_unclustered = argin.Results.show_unclustered;


% choose colormap
nclusters = max(cluster_index) - 1;
cmp = colormap(jet);
cmp = cmp(end:-1:1, :);
interval = floor(64/nclusters);

figure
if iscell(Orderofchannels) && numel(Orderofchannels) == 3
    clf;
    hold on
    for n = 1:nclusters
        color_map(n,:) = cmp(interval*n,:);
        scatter3(mydata(cluster_index == n,1), mydata(cluster_index == n,2),...
            mydata(cluster_index == n,3),10,color_map(n,:),'filled');
        cluster_name{n} = sprintf('cluster %d',n); % set legend
    end
    switch show_unclustered
        case 0
            % do not show unclustered cells
        case 1
            % show unclustered cells
            cluster_name{n+1} = 'unClustered cell';
            scatter3(mydata(cluster_index == n+1,1), mydata(cluster_index == n+1,2),...
                mydata(cluster_index == n+1,3),10,[208, 211, 212]/255,'filled');
        otherwise
            error('Invalid value for show_unclustered.')
    end
    legend(cluster_name);
    xlabel(Orderofchannels{1}); ylabel(Orderofchannels{2}); zlabel(Orderofchannels{3});
    hold off
    grid on
    savefig([filename '\visualize3D.fig'])
else
    error('The input Orderofchannels should be a cell array that contain 3 strings')
end           
