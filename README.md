# RainbowSeq
## Brief introduction
Core scripts for automatically gating in RainbowSeq.
The workflow for gating includes loading and pre-processing of image data (.tif) and FACS data (.fcs), clustering of pixels in image, mapping cells and pixels via "optimal transport", design for gating via "CART tree" and finally constructing the gate file (.xml) for cell sorting.

## Dependency
1. Matlab 2019a or 2019b
2. The supplied packages

## Usage (workflow of main function)
Here is an example. We provide biofilm image data (.tif), facs data (.fcs) , empty gate file (.xml) and main.m function for test.

### loading and pre-processing of image data (.tif)
Load image files (.tif) via `<imRead2>' and use `<iViewer>' to get the base fluorescene of biofilm to get the region of biofilm.  
```MATLAB
[~, colony] = colonyFL(img, base_fluorescnce)  
```
After get the region of biofilm. `<imgMat>` function could pre-process the image data and


