# RainbowSeq
## Brief introduction
Core scripts for automatically gating in RainbowSeq.
The workflow for gating includes loading and pre-processing of image data (.tif) and FACS data (.fcs), clustering of pixels in image, mapping cells and pixels via "optimal transport", design for gating via "CART tree" and finally constructing the gate file (.xml) for cell sorting.

## Dependency
1. Matlab 2019a or 2019b
2. The supplied packages （note: add them into the correct path)

## Usage 
Here is an example. We provide biofilm image data (.tif), facs data (.fcs) , empty gate file (.xml) and main.m function for test.
```MATLAB
matlab main.m
```  
Or open the main.m by matlab and run the pannels step by step.  

1. The main funtion (main.m) captures the pixels of biofilm from image data (.tif) and extract fluorescene of pixels from facs data by function `imgMat`. By this way, the fluorescene value of pixels and cells could be pre-processed by function `imgMat` and `facsMat` respectively.
2. Next pixels are clusterd into `nClusters` groups  via function`kmeans` according to the fluo. Hence, the cells from facs could be assigned to the pixels from image by `optimal transport (facsLabelAssign)` to obtain group information of corresponding pixels. Then function `Ttree` builds up a ***CART tree*** according to the group infromation and functions `splitpoint` and `saveTrees` parse the ***CART tree*** and transfrom the split methods of ***CART tree*** to gates in real facs space. Otherwise, function `calClusterInfo` calculates the metric (purity and yield) of each gate for further manually filtering. Finally, according to the final gates for each group of cells, `edixml2` edits the empty gate file (.xml) as template to generate the gating file used in cell sorting.

The successful runnning would generate the figures (.png or .fig),  gating strategy (.csv) and final gate file (.xml) shown in folder **.
