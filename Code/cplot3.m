function cplot3(mydata, Orderofchannels, cluster_index, filename)
%
% Display the clustered data in 3-D space.
% Usage: cplot3(mydata, Orderofchannels, cluster_index)
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
% calculated by function kmeans.]
%
% Written by Ping Shen @ Tsinghua University 
% Version 0.0.2. Created on Jul 15, 2019. Last modified on Jul 09, 2020.

% choose colormap
cmp = colormap(jet);
cmp = cmp(end:-1:1, :);
interval = floor(64/max(cluster_index));

% plot the 3D scatter

if iscell(Orderofchannels) && numel(Orderofchannels) == 3
else
    error('The input Orderofchannels should be a cell array that contain 3 strings')
end

figure
hold on
for n = 1:max(cluster_index)
    color_map(n,:) = cmp(interval*n,:);
    scatter3(mydata(cluster_index == n,1), mydata(cluster_index == n,2),mydata( cluster_index == n,3),10,color_map(n,:),'filled');
    cluster_name{n} = sprintf('cluster %d',n);
end
legend(cluster_name);
xlabel(Orderofchannels{1}); ylabel(Orderofchannels{2}); zlabel(Orderofchannels{3});
hold off
grid on
savefig([filename '\visualize3D.fig'])

clf;