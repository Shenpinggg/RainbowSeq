# RainbowSeq
## Brief introduction
Core scripts for automatically gating in RainbowSeq.
The workflow for gating includes loading and pre-processing of image data (.tif) and FACS data (.fcs), clustering of pixels in image, mapping cells and pixels via "optimal transport", design for gating via "CART tree" and finally constructing the gate file (.xml) for cell sorting.

## Dependency
1. Matlab 2019a or 2019b
2. The supplied packages （note: add them to search path of your matlab)

## Usage 
Here is an example. We provide biofilm image data (.tif), facs data (.fcs) , empty gate file (.xml) and main.m function.
```MATLAB
matlab main.m
```  
Or open the main.m by matlab and run the pannels step by step.  

1. The main funtion (main.m) captures the pixels of biofilm from image data (.tif) by set the threshold of fluorescence (`colonyFL`) and extract fluorescene of pixels from image  by function `imgMat`. <br>
***Correct biofilm detection (white line represents edge of biofilm)*** <br>
![image](https://github.com/Shenpinggg/RainbowSeq/blob/820b99da07e7ff9d3e09d5a91b7ade2a3e257e4e/Example/image_clustering/Test/colonyEdge.png)
***Poor detection (missing detection of biofilm interior)***<br>
![image](https://github.com/Shenpinggg/RainbowSeq/blob/820b99da07e7ff9d3e09d5a91b7ade2a3e257e4e/Example/image_clustering/Test/error_detection.png)
2. Pixels are clusterd into `nClusters` groups via function `kmeans` according to their fluorescence. Function `cmapping` could display the biofilm snapshot with clustering information for evaluation the result. <br> ***Biofilm clustered into 8 groups***<br>![image](https://github.com/Shenpinggg/RainbowSeq/blob/820b99da07e7ff9d3e09d5a91b7ade2a3e257e4e/Example/image_clustering/Test/mapping_to_biofilm.png) <br>
3. The fluorescence value of cells detected by FACS are extracted from facs file (.fcs) and pre-processed by function `facsMat`. Then the cells could be assigned to pixels from image by function `facsLabelAssign` (***optimal transport***) to get the grouping information of corresponding pixels. <br>
4. Then function `Ttree` builds up a ***CART tree*** according to the group infromation and functions `splitpoint` and `saveTrees2` parse the ***CART tree*** and transfrom the split methods of ***CART tree*** to gates in real facs space. Otherwise, function `calClusterInfo` calculates the metric (purity and yield) of each gate for further manually filtering. <br> ***Clusters' Abundance of pixels and cells*** <br> ***Valid method for gating*** <br> ![image]( https://github.com/Shenpinggg/RainbowSeq/blob/820b99da07e7ff9d3e09d5a91b7ade2a3e257e4e/Example/facs_clustering/Test/abundance.png)<br>***Poor method for gating***<br> ![image](https://github.com/Shenpinggg/RainbowSeq/blob/820b99da07e7ff9d3e09d5a91b7ade2a3e257e4e/Example/facs_clustering/Test/poor_result.png)<br>
5. Finally, according to the final gates for cell sorting, `edixml2` edits the empty gate file (.xml) as template to generate the gating file used in cell sorting. The final gate file (.xml) could be recognized for the BD FACSDiva software, directly loading all essential gates to the FACS machine.<br>

The successful runnning would generate the figures (.png or .fig),  gating strategy (.csv) and final gate file (.xml) shown in working directory `wd`.
