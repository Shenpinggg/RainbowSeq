# RainbowSeq
## Brief introduction
Core scripts for automatically gating in RainbowSeq.
The workflow for gating includes loading and pre-processing of image data (.tif) and FACS data (.fcs), clustering of biofilm pixels in image, mapping cells and pixels via "optimal transport", computing gating strategy via "CART tree" and finally constructing the gate file (.xml) for cell sorting.
For details about the algorithm, see out paper.

## Dependency
1. Matlab 2019a
2. The scripts in code directory （note: add these scripts to the PATH of your matlab)

## Usage 
Here is an example. We provide biofilm image data (.tif), facs data (.fcs) , empty gate file (.xml) and main.m function.
```MATLAB
matlab main.m (with default setting)
```  
Or open the main.m by matlab and run the pannels step by step (this is preferred, as several parameters need to be tuned manually).  

1. The main funtion (main.m) captures the pixels of biofilm from image data (.tif) by setting the threshold of fluorescence (`colonyFL`); next it extracts fluorescene of biofilm pixels from image by function `imgMat`. <br>
***Expected correct biofilm detection (white line represents edge of biofilm)*** <br>
![image](https://github.com/Shenpinggg/RainbowSeq/blob/820b99da07e7ff9d3e09d5a91b7ade2a3e257e4e/Example/image_clustering/Test/colonyEdge.png)
***Poor detection (missing detection of biofilm interior, need to set smaller threshold in this case)***<br>
![image](https://github.com/Shenpinggg/RainbowSeq/blob/820b99da07e7ff9d3e09d5a91b7ade2a3e257e4e/Example/image_clustering/Test/error_detection.png)
2. Pixels are clusterd into `nClusters` groups via function `kmeans` according to their fluorescence. Function `cmapping` could display the biofilm snapshot with group information to evaluate the performance of clustering. <br> ***One biofilm is perfectly clustered into 8 groups, giving rise to a rainbow structure***<br>![image](https://github.com/Shenpinggg/RainbowSeq/blob/820b99da07e7ff9d3e09d5a91b7ade2a3e257e4e/Example/image_clustering/Test/mapping_to_biofilm.png) <br>
3. The fluorescence value of cells detected by FACS are extracted from facs file (.fcs) and pre-processed by function `facsMat`. Then each cell could be assigned to one pixel from image by function `facsLabelAssign` (***optimal transport***); group label of the pixel was thus assigned to the corresponding cell. <br>
4. Then function `Ttree` builds up a ***CART tree*** according to the group infromation; next, functions `splitpoint` and `saveTrees2` parse the ***CART tree*** and convert the split boundaries of ***CART tree*** back to raw value (raw fluorescent value was normalized by function `facsMat`). Next, function `calClusterInfo` calculates the metric (purity and yield) of each gate. <br> 
***Clusters' Abundance of pixels (IMAGE) and cells (FACS)*** <br> 
***Good mapping (similar abundance between pixels and cells in the same cluster)*** <br> ![image](https://github.com/Shenpinggg/RainbowSeq/blob/820b99da07e7ff9d3e09d5a91b7ade2a3e257e4e/Example/facs_clustering/Test/abundance.png)<br>
***Poor mapping (quite different abundance between pixels and cells in the same cluster)***<br> ![image](https://github.com/Shenpinggg/RainbowSeq/blob/820b99da07e7ff9d3e09d5a91b7ade2a3e257e4e/Example/facs_clustering/Test/poor_result.png)<br>
5. Finally, according to the calculated gates, `edixml2` edits the empty gate file (.xml) as template to generate the gating file used in cell sorting. The final gate file (.xml) could be recognized for the BD FACSDiva software, directly loading all essential gates to the FACS machine.<br>
Note: the output gate file using this template file may be not recognized due to different versions of software or FACS machine; if so, prepare your own templaate file, do the same editing using `edixml2`.
The successful runnning would generate the figures (.png or .fig),  gating strategy (.csv) and final gate file (.xml) shown in working directory `wd`.
