%% loading image data
clear
wd = 'Example\';
cd(wd)
M_ph = imread('Test-PH.tif');
M_mCherry = imread('Test-mCherry.tif');
M_gfp = imread('Test-GFP.tif');
M_cfp = imread('Test-CFP.tif');
FACS_configure_path = 'FACS\configureFACS_biofilm.xlsx';
xmlFile = 'Template.xml';

imageRawMyData = {M_cfp, M_gfp, M_mCherry};
all_channels = {'CFP', 'GFP', 'mCherry'};
image_output_dir = 'image_clustering\Test\';
mkdir(image_output_dir);
rng default; % repeatable

%% Parameter to be set: number of clusters; channels of interest; thresholds of cells and pixels from FACS and image
nClusters = 8; % set cluster number
clustering_channels = {'CFP', 'GFP', 'mCherry'}; % channel for clustering
clustering_channelNum = [1,2,3];
clusterID = ['ClusterID', clustering_channels];

% set cutoff to reset outliers in image and cytometry
imgOutlierThr = [0.001, 0.001, 0.001];
facsOutlierThr = [0.001, 0.001, 0.001];
ncThr = [1.96, 1.96, 1.96];

%% identify pixel belonging to stained biofilm and calculate the distance to biofilm edge
line = [123.480932203389 2119.6281779661;3878 1675.50105932203]; % set the base line in the biofilm
colony = M_cfp >= 1200; % set the threshold of fluorescence to detect the pixel belonging to the biofilm
colony = medfilt2(colony,[5 5],'symmetric');
colony = imfill(colony,'holes');
[Dis, colonyE_v2] = distance2colonyEdge(colony, line, 'below'); % 'below' specifies the relative location of the abovementioned base line to biofilm in image
imshow(colonyE_v2) % check the performance of biofilm pixel identification, change the threshold of CFP if neccessary
saveas(gca,[image_output_dir 'colonyEdge.png'])

%% Image data pre-processing (extracting fluorescence and distance to edge of pixels)
[mydataNor, mydata] = imgMat(imageRawMyData, Dis, all_channels, colony,'outlier', imgOutlierThr);
distance = mydata(:,end);

%% K-means clustering of biofilm pixels into N groups according to fluorescent patterns
[idx_image, c_image] = kmeans(mydataNor(:, clustering_channelNum), nClusters, 'MaxIter',2000, 'Distance','cityblock');
idx_image = index_transition(idx_image, mydataNor(:, end));

%% Save the value of distance to biofilm edge for each cluster
for i = 1:nClusters
    cluster_mean_dis(i) = mean(mydata(idx_image == i,end))*-1;
end
cluster_mean_dis2edge = table(cluster_mean_dis, 'RowNames',{'Dis_Center2Edge'});
writetable(cluster_mean_dis2edge,[image_output_dir 'Cluster_mean_dis2edge.txt'],'Delimiter','tab','WriteRowNames',true);

%% Visualize clustering of biofilm pixels and the segmentation of biofilm
[coeff_np, score_np, latent_np, tsquared_np, explained_np, mu_np] = pca(mydataNor(:, clustering_channelNum),'Centered',false);
pca_plot(score_np, coeff_np, clustering_channels, image_output_dir, idx_image); % display clustering result of biofilm pixels via PCA
cmapping(idx_image, colony, image_output_dir); % biofilm segmentation according to the clustering result

%% FACS data pre-processing; mapping cells and pixels via "optimal transport"
facsMetaData = readtable(FACS_configure_path,'Sheet','Sheet1');

for dataNum = 1:height(facsMetaData)
    % get data
    biofilmFACSdata = facsMetaData.Biofilm{dataNum};
    ncFACSdata = facsMetaData.NC{dataNum};
    facs_output_dir = facsMetaData.Path{dataNum};
    mkdir(facs_output_dir);

    % process FACS data to extract information of each cell
    [myData_facs_Nor, myData_facs] = facsMat(biofilmFACSdata, ncFACSdata, ...
                                             all_channels, 'ncThr',ncThr, ...
                                             'outlierThr',facsOutlierThr);
                                             
    % remove all-negative cells (cells not stained by any one dye)
    for m = 1:length(clustering_channelNum)
        if m == 1
            nonBiofilmEntry = myData_facs_Nor(:, m) == 0;
        else
            nonBiofilmEntry = nonBiofilmEntry & (myData_facs_Nor(:, m) == 0);
        end
    end
    myData_facs_Nor = myData_facs_Nor(~nonBiofilmEntry, :);
    myData_facs = myData_facs(~nonBiofilmEntry, :);
    
    % Process the FACS data via Q normalization (optional)
    [qNorFACSdata,~] = qNorFACSimg(myData_facs_Nor(:,clustering_channelNum), mydataNor(:,clustering_channelNum), clustering_channels); 
    ForOTmapping_FACSdata = qNorFACSdata;
    
    % Prepare the image data for OT mapping
    dataImg = zeros(length(mydataNor), length(all_channels)+1);
    dataImg(:, 1:length(all_channels)) = mydataNor(:, 1:(end-1));
    dataImg(:, end) = idx_image;

    % OT mapping;
    % 'maxSubset' = 40 is to accelerate the computation;
    % 'OTsamplingFold' = 5 is to allow relaxation of the requirement that entry number between cells and pixels in each group should be all the same
    labeledFACS = facsLabelAssign(dataImg, ForOTmapping_FACSdata, ...
                                  'maxSubset',40, 'OTsamplingFold',5, ...
                                  'channels',clustering_channelNum);
    
    % Assign labeledFACS with pre-Qnorm data
    labeledFACS(:, clustering_channelNum) = myData_facs_Nor(:, clustering_channelNum);
    
end

% Visualize all cells in high-dimensional FL space
figure
cplot3(labeledFACS(:, clustering_channelNum), clustering_channels, labeledFACS(:, end), facs_output_dir);
close all

%% Visualize OT mapping result via PCA
[coeff_facs, score_facs, latent_facs, tsquared_facs, explained_facs, mu_facs] = pca(labeledFACS(:, clustering_channelNum),'Centered',false);
pca_plot(score_facs, coeff_facs, clustering_channels, facs_output_dir, labeledFACS(:, end));
close all

%% Visualize 10000 randomly sampled pixels in high-dimensional FL space
dataTmp = datasample(dataImg, 10000);
cplot3(dataTmp(:, clustering_channelNum), clustering_channels, dataTmp(:, end), image_output_dir);
close all

%% Save FACS data with the group information of each cell
myData_facs_index = [labeledFACS(:,end) myData_facs(:, clustering_channelNum)];
myData_facs_index = array2table(myData_facs_index, 'VariableNames', clusterID);
writetable(myData_facs_index, [wd '\facsRawData.csv'], 'Delimiter', 'tab')

%% Establish a CART tree for labeled FACS data (fluorescence as parameter, group information as label), based on which the FL space is splitted into cuboids
tree_facs = Ttree(labeledFACS(:,clustering_channelNum), labeledFACS(:,end), clustering_channels, nClusters+15, 'NumKfold', 10);

%% Use the CART tree model to assign each cell to relevant cuboid (cluster) in FL space
index_facs_cart = predict(tree_facs, labeledFACS(:,clustering_channelNum));
un = index_facs_cart ~= labeledFACS(:, end);
index_facs_cart_n = index_facs_cart;
index_facs_cart_n(un) = nClusters + 1; % label unclustered sample points
cplot3_2(labeledFACS(:,clustering_channelNum), clustering_channels, index_facs_cart_n, pwd,1);

%% Calculate the metric of each cuboid (cluster) in FL space
clusterInfo = calClusterInfo(labeledFACS(:,end), index_facs_cart, dataImg(:,end), facs_output_dir);

%% Convert the normalized data to raw fluorescence value in cytometry, thus we have the boundary information for each gate
[SplitPoints, boundTable] = splitpoint(tree_facs, myData_facs, myData_facs_Nor, clustering_channels);
[final_Decision, raw_tree, boundInfo, xmlCluster] = saveTree2(labeledFACS(:,end), myData_facs(:, clustering_channelNum), tree_facs, SplitPoints, ...
                                                              boundTable, clustering_channels, pwd, 'ViewTree', 1);

%% Edit the .xml template to generate the gating file used in final cell sorting
editxml2(xmlCluster, xmlFile, clustering_channels, wd);
