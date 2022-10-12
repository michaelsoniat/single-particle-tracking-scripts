# Tracking-Code

**Drift Correction**
To correct for sample drifting on the XY-stage plane, crop a quantum dot or other fluorphore that is fixed to the slide surface as a reference in Fiji. After saving as a tiff file, run IJF_Tracking.py by dragging it into Fiji and then hitting run which will output two files a results.tif and a .txt file with tracked results. If the program messes up on certain frames, the center can be manually added by going to Analyze>Tools>ROI Manager and adding an ROI to the list. The tracking program will automatically take into account this ROI as the center of the particle of interest. Following tracking, run registration-translatio-021115.ijm in Fiji. This will ask you to select the files of interest. First, pick the .txt file for the tracking of the background fluorophore. Then select the original movie which needs to be drift corrected. After run, this will output a new tiff file that has been drift corrected with the extention -registered which can be used for further data analysis. 

**Cropping**
The first step of tracking moving molecules is to crop out all of the molecules of interest in Fiji. Save the cropped images as tiff file so that all molecules can be tracked back to their position on the field. 

**Tracking**
After compiling a set of cropped molecules of interest, individually select and open those molecules in Fiji. Next, make sure that the frame scale is correct and crop to avoid any extraneous particles that might throw off the tracking program. Then, after saving that Analysis version of the file, run IJF_Tracking.py by dragging it into Fiji and then hitting run which will output two files a results.tif and a .txt file with tracked results. 

**Analysis**
Once several .txt files have been created from the tracking program, they should be analyzed with the ppfit_digest program in Matlab. Add the ppfit_digest.m file to the same folder as the .txt files along with the uigetfile2.m. In Matlab, in the prompt screen type ppfit_digest(# of kb/pixel, # of sec/frame). This will ask you to select the files of interest. Pick the .txt files by shift-clicking all of them. Then, it will ask you to select the regions of interest from each graph. Click and drag with the mouse to select the points of interest for each category and hit the letter key corresponding to that category to save those points as part of the category. The categories and their keys are below:

S = start plateau (the first frames in which molecule does not move). If molecule moves right away, simply select the first couple points to determine the start point.
D = decay slope (the regions where molecule moves in a straight line)
X = second slope (occasionally, molecules will pause then begin again, and this plots the slope of that movement)
F = final plateau (the final frames where molecule does not move (similar to S)

Finally, hit [A] to accept the data and save it or [Q] to disregard the plot and not save the data from that plot. 

This should create a set of .mat files. These can be called to Matlab individually and contain a structure called molec that contains all useful information such as molec.digest_rate and molec.x/molec.y.

Additionally, you can run the analyzeMats.m Matlab code to analyze all new .mat simultaneously and to create a new .mat file that combines all processivity and velocity values into a single file. To do this, run the analyzeMats.m file in Matlab. This will ask you to select the files of interest. Pick the .mat files by shift-clicking all of them. This should create two new set of .mat files which are  named Process and Rate.

**Please Cite:** Myler, L.R., Gallardo, I.F., Zhou, Y., Gong, F., Yang, S.H., Wold, M.S., Miller, K.M., Paull, T.T., and Finkelstein, I.J. (2016). Single-molecule imaging reveals
the mechanism of Exo1 regulation by single-stranded DNA binding proteins. Proc. Natl. Acad. Sci. USA 113, E1170–E1179.
