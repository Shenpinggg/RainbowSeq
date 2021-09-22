# RainbowSeq
## Brief introduction
Core scripts for automatically gating in RainbowSeq.
The workflow for gating includes loading and pre-processing of image data (.tif) and FACS data (.fcs), clustering of pixels in image, mapping cells and pixels via "optimal transport", design for gating via "CART tree" and finally constructing the gate file (.xml) for cell sorting.

## Dependency
1. Matlab 2019a or 2019b
2. The supplied packages

## Usage 
Here is an example. We provide biofilm image data (.tif), facs data (.fcs) , empty gate file (.xml) and main.m function for test. run
```MATLAB
matlab main.m
```  
Or open the main.m by matlab and run the pannels step by step.  

The main funtion (main.m) captures the pixels of biofilm from image data (.tif) and extract fluorescene of cells from facs data (.fcs). By this way, the fluorescene value of pixels and cells could be pre-processed by function `imgMat` and `facsMat` respectively.  
Next pixels are clusterd into N groups via function`kmeans`. Then the cells from facs could be assigned to the pixels from image by `optimal transport (facsLabelAssign)` to obtain group information of corresponding pixels. Hence, the function `Ttree` builds up a ***CART tree*** according to the group infromation of cells for gate design and functions splitpoint. Finally, according to the final gate for each group of cells, `edixml2` edits the empty gate file (.xml) as template to generate the gating file used in cell sorting
