# RainbowSeq
## brief introduction
Core scripts for automatically gating in RainbowSeq.
The workflow for gating includes loading and pre-processing of image data (.tif) and FACS data (.fcs), clustering of pixels in image, mapping cells and pixels via "optimal transport", design for gating via "CART tree" and finally constructing the gate file (.xml) for cell sorting.

## Dependency
1. Matlab 2019a or 2019b
2. The supplied packages

## Usage (workflow of main function)
Here is an example. We provide biofilm image files (.tif), facs data (.fcs) , empty gate file (.xml) and main.m function for test.
First load the image files and detect the region of biofilm.
```MATLAB
ph = imRead2('img1.tif')

