# RainbowSeq
## Brief introduction
MATLAB scripts for automatically compute gating strategy in RainbowSeq.
The workflow for gating includes loading and pre-processing of image data (.tif) and FACS data (.fcs), clustering of biofilm pixels in image, mapping cells and pixels via "optimal transport", computing gating strategy via "CART tree" and finally constructing the gate file (.xml) for cell sorting.
For details about the algorithm, see our paper.

## Dependency
1. MATLAB 2019a
2. The scripts in Code directory （note: add these scripts to the PATH of your MATLAB)

## Input files (in current folder) <br>
Image files: `/Example/Test-PH/CFP/GFP/mCherry.tif` <br>
FACS files: `/Example/FACS/NC.fcs` (negative control) & `/Example/FACS/biofilm.fcs` (experiment data) <br>
Configure of FACS file: `/Example/FACS/configureFACS_biofilm.xlsx` <br>
![image](https://github.com/Shenpinggg/RainbowSeq/blob/86ff87142df7b3b3be7f856fbe63d99190bbf8b6/Example/FACS/configure_FAS_file.png) <br>
Biofilm: Path of biofilm facs file (biofilm.fcs in example) <br>
NC: Path of negative control facs file (NC.fcs in example) <br>
Path: Path stored with figures involved in FACS data overview, FACS calibration and mapping <br>
Template xml file: `/Example/Template.xml` export from FACS machine <br>
***Path of above input files could be modified in 1st section of main.m***

## Usage 
Here is an example. We provide biofilm image data (.tif), facs data (.fcs) , empty gate file (.xml) and main.m function.
```MATLAB
matlab main.m (with default setting)
```  
Or open the main.m by matlab and run the pannels step by step (this is preferred, as several parameters need to be customized).  

1. The main funtion (main.m) captures the pixels of biofilm from image data (.tif) by setting the threshold of fluorescence (`colonyFL`); next, it extracts fluorescene of biofilm pixels from image by function `imgMat`. <br>
***Expected correct biofilm detection (white line represents edge of biofilm)*** <br>
![image](https://github.com/Shenpinggg/RainbowSeq/blob/820b99da07e7ff9d3e09d5a91b7ade2a3e257e4e/Example/image_clustering/Test/colonyEdge.png)
***Poor detection (fail to detect biofilm interior, need to set smaller threshold in this case)***<br>
![image](https://github.com/Shenpinggg/RainbowSeq/blob/820b99da07e7ff9d3e09d5a91b7ade2a3e257e4e/Example/image_clustering/Test/error_detection.png)
2. Pixels are clusterd into `nClusters` groups via function `kmeans` according to their fluorescence. Function `cmapping` could display the biofilm snapshot with group information to evaluate the performance of clustering. <br> ***One biofilm is perfectly clustered into 8 groups, giving rise to a rainbow structure***<br>![image](https://github.com/Shenpinggg/RainbowSeq/blob/820b99da07e7ff9d3e09d5a91b7ade2a3e257e4e/Example/image_clustering/Test/mapping_to_biofilm.png) <br>
3. The fluorescence value of cells detected by FACS are extracted from facs file (.fcs) and pre-processed by function `facsMat`. Then each cell could be assigned to one pixel from image by function `facsLabelAssign` (***optimal transport***); group label of the pixel was thus assigned to the corresponding cell. <br>
4. Then function `Ttree` builds up a ***CART tree*** according to the group labels (label) and fluorescent values (feature) of all cells; next, functions `splitpoint` and `saveTrees2` parse the ***CART tree*** and convert the split boundaries of ***CART tree*** back to raw fluorescent value (raw fluorescent value was normalized by function `facsMat` at previous step). Next, function `calClusterInfo` calculates the metric (purity and yield) of each gate. <br> 
***Clusters' Abundance of pixels (IMAGE) and cells (FACS)*** <br> 
Note that only focus on abundance of pixels (blue bar below) and cells (orange bar below), which suggests the quality of mapping.
The yellow bar below specifies the final abundance of each cluster within the relevant gates from CART-tree. This metric is sometimes different from the relevant pixel or cell abundance, as the geometry of cells of such cluster in high-dimensional fluorescent space cannot be completely captured by the rectangle gates derived from CART-tree. This does not interfere with the accuracy of the method.
***Good mapping (similar abundance between pixels and cells in the same cluster)*** <br> ![image](https://github.com/Shenpinggg/RainbowSeq/blob/820b99da07e7ff9d3e09d5a91b7ade2a3e257e4e/Example/facs_clustering/Test/abundance.png)<br>
***Poor mapping (quite different abundance between pixels and cells in the same cluster)***<br> ![image](https://github.com/Shenpinggg/RainbowSeq/blob/820b99da07e7ff9d3e09d5a91b7ade2a3e257e4e/Example/facs_clustering/Test/poor_result.png)<br>
5. Finally, according to the calculated gates, `edixml2` edits the empty gate file (.xml) as template to generate the gating file used in cell sorting. The final gate file (.xml) could be recognized by the BD FACSDiva software, directly loading all essential gates to the FACS machine.<br>
Note: the output gate file using this template file may be not recognized due to different versions of software or FACS machine; if so, prepare your own templaate file, do the same editing using `edixml2`.<br>
## Output
***Image clustering*** <br>
Path of folder: `/Example/image_clustering/Test/`
1. Clustering mapped to biofilm `mapping_to_biofilm.png` <br>
2. Distance (micrometer) of clusters to biofilm edge `Cluster_mean_dis2edge.txt`<br>
3. View fluorescence of clustered pixels in 3D space `visualize3D.png` <br>

***Map cells to clustered pixels*** <br>
Path of folder: `/Example/facs_clustering/Test/`
1. Fluoresence of cells of different channels detected by FACS `RainbowSeq/Example/facsRawData.csv`
2. View fluorescence of mapped cells in 3D space  `visualize3D.png` (The topology of fluorescence from cells and pixels are similar in good mapping) <br>
3. Abundance between pixels and cells in each cluster `abundance.png`<br> 

***Sorting strategy***
Path of folder: `/Example/`
1. Split strategy of ***CART tree*** (rawdata of the tree) to generate the below gating strategy files `treeDecision.csv`
2. To facilitate the usage of gate file, the program also generates a sereis of .xlsx files to specify the relevant gates for each cluster.
The .xlsx file is named after custerX_final_gating_strategy.csv, such files will be generated according to the number of groups during k-means clustering of biofilm pixels.
The screenshot of `RainbowSeq/Example/Cluster1_final_gating_strategy.csv` <br>
![image](https://github.com/Shenpinggg/RainbowSeq/blob/bbb5695ff8bcc4ac094641babbec62a0ac114cc4/Example/final_gating_strategy_Cluster1.png)<br>
3. These .xlsx files are essential to guide the usage of the final gate file (.xml) `/Example/calculated_gate.xml` for cell sorting.
In each file, we have a series of gates targeting one group of cells. As shown above, each row specify such a (combined) gate.
Purity suggests the abundance of cells belonging to the desired group among all cells in this gate; yield suggests the abundance of desired cells in this gate among all cells; relative yield suggests the abundance of desired cells in this gate among all cells belonging to this group.
In some cases, we have gate1 and gate2 in one row. This structure suggests that combining gate1 and gate2 will generate one proper gate for relevant cells. In fact, we have editted the gate file (.xml) to set gate2 as the son gate of gate1 in such cases. Hence, just choose gate2 in FACS software in such cases; if only gate1 is available, choose gate1. For example, we can choose gateX and gateY to sort the group of cells shown in this example.
Another tip is to manually filter some gates with low purity or yield. This improves both the accuracy of the method and saves time. Usually, gates with low purity (< 80%) and low yield (relative yield < 1%) were excluded.
