# ibb_toolkit
Iowa Brain-Behavior Toolkit: An open-source MATLAB tool for inferential and predictive modeling of imaging-behavior and lesion-deficit relationships

The /tutorials folder contains MATLAB Live notebooks illustrating how to perform different kinds of modeling analyses using the toolkit. These tutorials illustrate how to import different kinds of datasets, how to specify and run different modeling analyses, and how to output and evaluate analysis results. 

The /manual folder contains a User Manual that provides an overview of toolkit functionality. 

To invoke the GUI, navigate to the /ibb_toolkit directory in MATLAB, and type run_modeling_gui into the command window. 

NOTE - Running PLSR/PLS-DA without stratification will break in R2023b (and possibly other versions after R2022b) due to a bug in plsregress() that prevents using both the 'CV' and 'MCReps' name-value pair arguments with integer values. I have reported the bug to MATLAB and am awaiting resolution. Please either (a) use MATLAB R2022b for PLSR/PLS-DA analyses if stratification is not wanted, or (b) use stratification for PLSR/PLS-DA analyses. 
