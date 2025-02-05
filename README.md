# ibb_toolkit
Iowa Brain-Behavior Toolkit: An open-source MATLAB tool for inferential and predictive modeling of imaging-behavior and lesion-deficit relationships

The /tutorials folder contains MATLAB Live notebooks illustrating how to perform different kinds of modeling analyses using the toolkit. These tutorials illustrate how to import different kinds of datasets, how to specify and run different modeling analyses, and how to output and evaluate analysis results. 

The /manual folder contains a User Manual that provides an overview of toolkit functionality. 

To invoke the GUI, navigate to the /ibb_toolkit directory in MATLAB, and type run_modeling_gui into the command window. 

NOTE - Permutation testing for mass-univariate GLMs that can utilize covariates (ordinary least squares, logistic regression multinomial regression) can be relatively slow, but parametric p-values can also be obtained from these models without permutation testing; however, this precludes use of the cFWE procedure and requires correction across all voxels, which may be more conservative. 

NOTE - Running PLSR/PLS-DA without stratification will break in R2023b due to a bug in plsregress() that prevents using both the 'CV' and 'MCReps' name-value pair arguments with integer values. I have reported the bug to MATLAB, and they basically said "it's fixed in later versions, so just use those". Please either (a) use MATLAB R2022b or a version later than R2023b for PLSR/PLS-DA analyses if stratification is not wanted, or (b) use stratification for PLSR/PLS-DA analyses if you are running in R2023b. 
