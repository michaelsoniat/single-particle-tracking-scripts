This is a collection of legacy scripts for rudimentary particle tracking and trajectory analysis with the [FIJI](https://fiji.sc/) image processing package. Particles are tracked in one dimension by fitting the point spread function to a 2D Gaussian. The resulting trajectories are analyzed via a series of MATLAB scripts, as described below.

Please see [Alternative Packages](#alternative-packages) for newer and more general particle tracking workflows.

<!-- markdown-toc start - Don't edit this section. Run M-x markdown-toc-refresh-toc -->
**Table of Contents**

- [Table of Contents](#table-of-contents)
- [Using this package](#using-this-package)
    - [Drift Correction](#drift-correction)
    - [Cropping](#cropping)
    - [Tracking](#tracking)
    - [Analysis](#analysis)
- [Citations](#citations)
- [Alternative Packages](#alternative-packages)

<!-- markdown-toc end -->


# Using this package #

The workflow below captures how we analyze the motion of fluorescent proteins on laminar flow-stretched DNA molecules.

## Drift Correction ##

For long-term imaging experiments, the microscope stage can drift significantly in the x-y plane. To correct for sample drifting, select several fluorophores that are fixed to the flowcell surface. Crop these reference particles as a tiff stack in FIJI and run `IJF_Tracking.py` by dragging it into FIJI and then hitting "run." This script outputs two files: (1) results.tif and (2) results.txt, which has x-y coordinates for the particle position in each TIFF frame. 

Note: If the script loses the particle position on some frames, the center can be manually added by going to Analyze>Tools>ROI Manager and adding a point ROI to the list. The tracking script will automatically fit a 2D Gaussian near this ROI for that frame.

After obtaining x-y coordinates for the stationary particle, run `registration-translatio-021115.ijm` in FIJI. This will ask you to select the files of interest. First, pick the txt file containing the x-y coordinates for each frame of the stationary fluorophore. Then select the TIFF stack that needs to be drift-corrected. This script will output a new tiff file that has been drift-corrected for further data analysis. 

## Cropping ##

The goal of this step is to reduce the chance that the particle tracking script will lose its focus on the particle of interest and start tracking another particle that may be adjacent to it. 

Crop out the molecules of interest in FIJI. Make sure that the frame scale is correct and that the ROI avoids all particles that might throw off the tracking program. Save the cropped images as separate TIFF stacks for the next step.

## Tracking ##

Select and open the cropped TIFF stacks in FIJI. Run `IJF_Tracking.py` by dragging it into FIJI and then hitting run which will output two files, a results.tif and a txt file with the x-y coordinates for each TIFF frame. 

## Analysis ##

We use [MATLAB](https://www.mathworks.com/products/matlab.html) scripts with a rudimentary user interface to extract useful information from particle trajectories.

The `ppfit_digest.m` script fits a subset of each particle trajectory to a piecewise polynomial. This script also uses [uigetfile2](https://www.mathworks.com/matlabcentral/fileexchange/9254-uigetfile2) to improve the file selection user interface. I suggest you copy `uigetfile2.m` into the same folder as `ppfit_digest.m`. 

In MATLAB, run `ppfit_digest(kb_per_pixel, sec_per_frame)`. The first parameter is a conversion factor between kilobases and pixels and the second converts between seconds and frames. These settings depend on the microscope and image acquisition settings.

When prompted, select the TXT files with the drift-corrected trajectories. The script will then display the particle trajectory on a graph. Click and drag with the mouse to select the points of interest that define a distinct sub-trajectory (i.e., a paused state or a translocating state). The script defines the following key shortcuts for fitting these trajectories to piecewise polynomials:

``` 
s = start (the first frames in which the molecule does not move). If molecule moves right away, simply select the first couple points to determine the starting point.

d = downward slope (the regions where molecule moves in an approximately straight line, e.g., DNA unwinding by a helicase)
x = second slope (occasionally, molecules will pause then begin again, and this plots the slope of that movement)
f = final plateau (the final frames where molecule does not move (similar to [s])

a = accept the data and save it
q = disregard this trajectory and move on to the next one

```
Fitting trajectories to polynomials will create a set of MAT files that contain a structure that contains the piecewise polynomial fits.

Additionally, you can run `analyzeMats.m` to analyze all new MAT files simultaneously. Running this script will ask you to select the files of interest. Pick the MAT files by shift-clicking all of them. The output is a MAT file with summary statistics about the analyzed trajectories.

# Citations #

If you find this script useful, please cite:

Finkelstein IJ, Visnapuu ML, Greene EC. Single-molecule imaging reveals mechanisms of protein disruption by a DNA translocase. Nature. 2010 Dec 16;468(7326):983-7. doi: [10.1038/nature09561](https://doi.org/10.1038/nature09561). Epub 2010 Nov 24. PMID: 21107319; PMCID: [PMC3230117](https://pubmed.ncbi.nlm.nih.gov/21107319/).

# Alternative Packages #
This script was developed to track single-particle trajectories in one dimension and has been superseded by much more powerful packages. The incomplete list below summarizes a few resources to get the user started.

**Mars-Molecule Archive Suite**

Citation: Nadia M Huisjes Thomas M Retzer Matthias J Scherr Rohit Agarwal Lional Rajappa Barbara Safaric Anita Minnen Karl E Duderstadt (2022) Mars, a molecule archive suite for reproducible analysis and reporting of single-molecule properties from bioimages. eLife, 11:e75899. [doi: 10.7554/eLife.75899](https://doi.org/10.7554/eLife.75899)

**Trackmate**

Citation: Tinevez, J.-Y., Perry, N., Schindelin, J., Hoopes, G. M., Reynolds, G. D., Laplantine, E., … Eliceiri, K. W. (2017). TrackMate: An open and extensible platform for single-particle tracking. Methods, 115, 80–90. doi:10.1016/j.ymeth.2016.09.016
