# How to work with Modified GALFIT Python3 wrapper

Author: Bingcheng Jin

Date: 2024.1.1

## GALFITM Introduction

GALFIT is a two-dimensional (2D) fitting algorithm developed by Peng et al. (2002, 2010), which uses the Levenberg-Marquardt algorithm to minimize the $\chi^2$ between a galaxy image and the PSF-convolved model image.

GALFIT is available at https://users.obs.carnegiescience.edu/peng/work/galfit/galfit.html.

GALFITM is a modified version of GALFIT, which can simultaneously fit multiple wavelength images to obtain wavelength-dependent morphology properties. It performs well especially for faint galaxies with low surface brightness.

GALFITM is available at https://www.nottingham.ac.uk/astronomy/megamorph/.

## Installation and Set-up

To fully reproduce our fit, please install the latest version of GALFITM v1.4.4,
which enable the use of higher-order functions (bending, fourier modes, rotation) from GALFIT3.

To run GALFITM in command line, we recommand you alias the software downloaded
from https://www.nottingham.ac.uk/astronomy/megamorph/.

For example, in Bash:
```
alias galfitm='~/Downloads/galfitm-1.4.4-osx'
```


### Python Environment

```
conda create -n <env_name> python=3.9
conda activate <env_name>
```

### Additional packages and software

Latest version of Astropy and Photutils is needed.
```
pip install astropy==5.3.1
pip install photutils==1.9.0
```


## Stage 1: Source Detection ([Notebook](1catsrc.ipynb))

### [Merge Image](utils/merge.py)
Merge 4 LW bands SCI_BKSUB extension to produce an stacked image for source detection. Their weight are determined by the error map (inverse of the ERR extension).

### Perform Convolution
After careful consideration, the convolution kernel are chosen to be little bigger than the FWHM of the NIRCam LW PSFs.

### Source Detection
In order to detect as many galaxies as possible, the threshold for source detection is set around 1.05 $\sigma$.

### Crossmatch with Catalog (Stefanon et al. 2017)
The matching results are not perfect as expected, which requires manual confirmation:
* [Manual Matching](manualcheck.py)
```
python manualcheck.py -i <img_name>
```
or simply run the code, it will ask you which image to process.
* [Implementary Notebook](manualcomplement.ipynb)

## Stage 2: Photometry Preperation ([Notebook](2prephot.ipynb))

### Create Source Catalog using photutils
This step provides initial guess for following GALFIT process.

### Generate pregalfit.fits for following steps
The initial parameters, segmentation map and image list will be stored in the file for convenience.

## Stage 3: Cutout Images 

### Create [SegmentCut Class](utils/segmentcut.py) Object
SegmentCut is designed for convenience of cutting images with input cutout size, centroid position, targeted labels etc.

### Image Cutout
To alleviate the stress of dealing with large numbers of data, cuting out images is necessary. The size will be based on the area segmented by source detection process.
* Ensemble [Image Cutting](cutout.py)
```
python cutout.py -i <img_name> -o
```
or simply run the code, it will ask you which image to process.

### Free to recut 
* [Manual Cutting](3recut.ipynb)

## Stage 4: Run GalfitM ([Notebook](4dothefit.ipynb))

### Create [Galfit Class](galfitclass.py) Object
Galfit Class is designed for better experience on the Python Interface. It will help manipulating the galfit configuration (including the overall setup and components).

### Generate Feedmes, Link Contraints
This is also covered by the Galfit Class.

### Run GALFITM with Command Lines
Command lines are also embedded in the wrapper, providing better user experience in coding.

* [Ensemble Fitting Python Code](runfit.py)
```
python runfit.py -i <img_name> -o
```
or simply run the code, it will ask you which image to process.

### Stage 5: Post-pipeline Analysis ([Notebook](collect.ipynb))
collect output files from all images and export to [ecsv file](result.ecsv)
