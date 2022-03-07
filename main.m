%% loading image data
clear
wd = 'Example\';
cd(wd)
M_ph = imread('Test-PH.tif');
M_mCherry = imread('Test-mCherry.tif');
M_gfp = imread('Test-GFP.tif');
M_cfp = imread('Test-CFP.tif');

imageRawMyData = {M_cfp, M_gfp, M_mCherry};
all_channels = {'CFP', 'GFP', 'mCherry'};
image_output_dir = 'image_clustering\Test\';
mkdir(image_output_dir);
rng default; % repeatable
%% Parameter: number of clusters; channels of interest; thresholds of cells and pixels from FACS and image
nClusters = 8; % set cluster number
clustering_channels = {'CFP', 'GFP', 'mCherry'}; % channel for clustering
clustering_channelNum = [1,2,3];
clusterID = ['ClusterID', clustering_channels];
% remove outliers in image
imgOutlierThr = [0.001, 0.001, 0.001];
facsOutlierThr = [0.01, 0.01, 0.0001];
ncThr = [1.96, 1.96, 1.96];
%% get the colony and the distance to the edge
line = [123.480932203389 2119.6281779661;3878 1675.50105932203]; % set the base line in biofilm
colony = M_cfp >= 1200; % set the threshold of fluorescence to detect the biofilm
colony = medfilt2(colony,[5 5],'symmetric'); 
colony = imfill(colony,'holes');
[Dis, colonyE_v2] = distance2colonyEdge(colony, line, 'below');
imshow(colonyE_v2)
saveas(gca,[image_output_dir 'colonyEdge.png'])
%% Image data pre-processing (extracting fluorescence and distance to edge of pixels)
[mydataNor, mydata] = imgMat(imageRawMyData, Dis, all_channels, colony,'outlier', imgOutlierThr, 'position_coeff', 2);
distance = mydata(:,end);
%% K-means clustering of pixels in 3D fluorescent space into N groups
[idx_image,c_image] = kmeans(mydataNor(:, clustering_channelNum), nClusters, 'MaxIter',2000, 'Distance','cityblock');
idx_image = index_transition(idx_image, mydataNor(:, end));
%% save the distance of clusters 
for i = 1:nClusters
    cluster_mean_dis(i) = mean(mydata(idx_image == i,end))*-1;
end
cluster_mean_dis2edge = table(cluster_mean_dis, 'RowNames',{'Dis_Center2Edge'});
writetable(cluster_mean_dis2edge,[image_output_dir 'Cluster_mean_dis2edge.txt'],'Delimiter','tab','WriteRowNames',true);
%% Visualize image clustering in FL space 
[coeff_np, score_np, latent_np, tsquared_np, explained_np, mu_np] = pca(mydataNor(:, clustering_channelNum),'Centered',false);
pca_plot(score_np, coeff_np, clustering_channels, image_output_dir, idx_image); % display in PCA space
cmapping(idx_image, colony, image_output_dir); % map back to the biofilm
%% FACS data pre-processing; Mapping cells and pixels via "optimal transport"
facsMetaData = readtable('FACS\configureFACS_biofilm.xlsx','Sheet','Sheet1');

for dataNum = 1:height(facsMetaData)
    % get data
    biofilmFACSdata = facsMetaData.Biofilm{dataNum};
    ncFACSdata = facsMetaData.NC{dataNum};
    facs_output_dir = facsMetaData.Path{dataNum};
    mkdir(facs_output_dir);

    % process data
    [myData_facs_Nor, myData_facs] = facsMat(biofilmFACSdata, ncFACSdata, ...
                                             all_channels, 'ncThr',ncThr, ...
                                             'outlierThr',facsOutlierThr);    
    % remove all-negative cells, possibly contaminated
     for m = 1:length(clustering_channelNum)
         if m == 1
             nonBiofilmEntry = myData_facs_Nor(:, m) == 0;
         else
             nonBiofilmEntry = nonBiofilmEntry & (myData_facs_Nor(:, m) == 0);
         end
     end
    myData_facs_Nor = myData_facs_Nor(~nonBiofilmEntry, :);
    myData_facs = myData_facs(~nonBiofilmEntry, :);
    
    % Process the FACS data via Q normalization
    [qNorFACSdata,~] = qNorFACSimg(myData_facs_Nor(:,clustering_channelNum), mydataNor(:,clustering_channelNum), clustering_channels); 
    ForOTmapping_FACSdata = qNorFACSdata;
    
    % Prepare the image data for OT mapping
    dataImg = zeros(length(mydataNor), length(all_channels)+1);
    dataImg(:, 1:length(all_channels)) = mydataNor(:, 1:(end-1));
    dataImg(:, end) = idx_image;
    
    % Present the comparison of distribution between FACS and image for each channel
    img2facsSim(mydataNor(:, 1:(end-1)), myData_facs_Nor, ...
                all_channels, facs_output_dir);

    % OT mapping
    labeledFACS = facsLabelAssign(dataImg, ForOTmapping_FACSdata, ...
                                  'maxSubset',40, 'OTsamplingFold',5,...
                                  'channels',clustering_channelNum);
    
    % Assign labeledFACS with pre-Qnorm data
    labeledFACS(:, clustering_channelNum) = myData_facs_Nor(:, clustering_channelNum);
    
    % result of FACS
    groupLabels = unique(labeledFACS(:, end));
    cmp = colormap(jet); cmp = cmp(end:-1:1, :); interval = floor(64/max(groupLabels));
    conditions = {};
    binNum = 50;
    hold on
    for group = 1:max(groupLabels) % each cluster
        color_map(group,:) = cmp(interval*group,:);
        thisGroup = labeledFACS(labeledFACS(:,end) == group, end-1);
        histogram(thisGroup, binNum, 'EdgeColor',color_map(group,:), ...
                  'DisplayStyle', 'stairs', 'LineWidth',1, 'Normalization','probability');
        conditions{end+1} = ['cluster', num2str(group)];
    end
    hold off
    alpha(0.6);
    lgd = legend(conditions, 'Location','northeast','FontSize',10);
    xlabel('mCherry of cells','Fontsize',18)
    ylabel('Probability','Fontsize',18)
    filename = [facs_output_dir, '\groupsByHistFACS.png'];
    saveas(gcf, filename, 'png')
    
    % result of image
    conditions = {};
    clf();
    hold on
    for group = 1:max(groupLabels) % each cluster
        color_map(group,:) = cmp(interval*group,:);
        thisGroup = dataImg(dataImg(:,end) == group, end-1);
        histogram(thisGroup, binNum, 'EdgeColor',color_map(group,:), ...
                  'DisplayStyle', 'stairs', 'LineWidth',1, 'Normalization','probability');
        conditions{end+1} = ['cluster', num2str(group)];
    end
    hold off
    alpha(0.6);
    lgd = legend(conditions, 'Location','northeast','FontSize',10);
    xlabel('mCherry of pixels','Fontsize',18)
    ylabel('Probability','Fontsize',18)
    filename = [facs_output_dir, '\groupsByHistImg.png'];
    saveas(gcf, filename, 'png')
end
% Visualize FACS in FL space
figure
cplot3(labeledFACS(:, clustering_channelNum), clustering_channels, labeledFACS(:, end), facs_output_dir);
close all

%% Visualize mapping of FACS in PCA space
[coeff_facs, score_facs, latent_facs, tsquared_facs, explained_facs, mu_facs] = pca(labeledFACS(:, clustering_channelNum),'Centered',false);
pca_plot(score_facs, coeff_facs, clustering_channels, facs_output_dir, labeledFACS(:, end));
close all
%% Visualize sampled image in FL space
dataTmp = datasample(dataImg, 10000);
cplot3(dataTmp(:, clustering_channelNum), clustering_channels, dataTmp(:, end), image_output_dir);
close all
%% save the facs data as csv file 
myData_facs_index = [labeledFACS(:,end) myData_facs(:, clustering_channelNum)];
myData_facs_index = array2table(myData_facs_index, 'VariableNames', clusterID);
writetable(myData_facs_index, [wd '\facsRawData.csv'], 'Delimiter', 'tab')
%% establish a CART tree of FACS data
tree_facs = Ttree(labeledFACS(:,clustering_channelNum), labeledFACS(:,end), clustering_channels, nClusters+15, 'NumKfold', 10);
%% use the CART tree model to assign each entry to relevant cluster
index_facs_cart = predict(tree_facs, labeledFACS(:,clustering_channelNum));
un = index_facs_cart ~= labeledFACS(:, end);
index_facs_cart_n = index_facs_cart;
index_facs_cart_n(un) = nClusters + 1; % label unclustered sample points
cplot3_2(labeledFACS(:,clustering_channelNum), clustering_channels, index_facs_cart_n, pwd,1);
%% calculate the metric of each cluster
clusterInfo = calClusterInfo(labeledFACS(:,end), index_facs_cart, dataImg(:,end), facs_output_dir);
%% transform the normalized data to real FACS space and get the split method
[SplitPoints, boundTable] = splitpoint(tree_facs, myData_facs, myData_facs_Nor, clustering_channels);
[final_Decision, raw_tree, boundInfo, xmlCluster] = saveTree2(labeledFACS(:,end), myData_facs(:, clustering_channelNum), tree_facs, SplitPoints, boundTable, clustering_channels, pwd, 'ViewTree', 1);
%% edit the xml template to generate the gating file used in final cell sorting
xmlFile = 'Template.xml';
editxml2(xmlCluster, xmlFile, clustering_channels, wd);