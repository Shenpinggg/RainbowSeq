function img2facsSim(imgData, facsData, channelsOfInterest, filename)

% This function is to compare the histogram of entries (pixels and cells)
% belonging to image and FACS.
%
% Usage: img2facsSim(myDataNor, myData_facs_Nor, channelsOfInterest, filename)
%
% The input IMGDATA should be a n*p array, where p is number of channels
% and n is number of pixels.
%
% The input FACSDATA should be a m*p array, where p is number of channels
% and m is number of cells.
%
% The input channelsOfInterest is a cell array, specifying the channels of
% interest. Accepted channels include 'CFP', 'GFP', 'mCherry' and 'APC'. It
% should have p channels.
%
% The input filename is the integrate path as well as the prefix of the
% output plot, such as 'E:\Data\D0500\' (no prefix) or
% 'E:\Data\D0500\lowVoltage'
%
% Output
% The program will save p plots in given path. Each of them is two
% histogram of both image and FACS for a particular channel. The Pearson
% correlation coefficient is also presented in the plot.
%
% Written by Tianmin Wang @ Tsinghua University
% Version 0.0.1. Created on Aug 11, 2019. Last modified on Aug 12, 2019.

% check data dimension consistency
if ndims(imgData) ~= 2 || ndims(facsData) ~= 2
    error('Input dataset should be 2D array!')
end

% check channel number consistency
[Ny_img, Nx_img] = size(imgData);
[Ny_facs, Nx_facs] = size(facsData);

if Nx_img ~= numel(channelsOfInterest) || Nx_facs ~= numel(channelsOfInterest)
    error('Number of channels is not consistent!')
end

% get imgae and facs datasets with the same amount of entries
datasetSize = min(Ny_img, Ny_facs);
myDataNor_sample = datasample(imgData, datasetSize);
myData_facs_Nor_sample = datasample(facsData, datasetSize);

% color
colors = struct;
colors.CFP = 'b';
colors.GFP = 'g';
colors.mCherry = 'm';
colors.APC = 'r';
edges = [0:0.02:1]; % histogram boundary

% comparison between image and facs for each channel
for m = 1:numel(channelsOfInterest)
    clf;
    channel = channelsOfInterest{m};
    color = colors.(channel);
    channelmIMG = sort(myDataNor_sample(:, m),'descend');
    channelmFACS = sort(myData_facs_Nor_sample(:, m),'descend');
    R = corrcoef(channelmIMG, channelmFACS); % correlation
    corrP = R(1,2);
    hold on
    histogram(channelmIMG, edges, ...
                  'DisplayStyle','stairs', 'EdgeColor','k');
    histogram(channelmFACS, edges, ...
                  'DisplayStyle','stairs', 'EdgeColor',color);
    hold off
    alpha(0.3);
    title(['Channel: ', channel, '; PCC: ', num2str(corrP)]);
    xlabel('Fluorescence', 'Fontsize',18)
    ylabel('Number of entries', 'Fontsize',18)
    lgd = legend({' Image: pixels', ' FACS: cells'}, 'Location','northeast');
    lgd.FontSize = 14;
    saveas(gcf, [filename, 'channel_', channel, '.png'])
end

clf;