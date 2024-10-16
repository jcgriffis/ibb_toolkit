README

***Notes about graphical user interface (GUI)***

To run the GUI, type run_modeling_gui into the MATLAB command window while in the main toolkit directory or after adding it to your path.

To facilitate file and directory selection, you can create a file named "default_paths.mat" in this folder that contains the following variables:

 - csv_dir - character array or string containing path to directory cpntaining the CSV file(s) with behavioral data
 - data_dir - character array or string containing path to directory containing patient predictor files (e.g., lesion mask NIFTI images)
 - mask_file - character array or string containing path to the brain mask NIFTI file (for analyses using voxel-wise data)
 - parcel_file - character array or string containing path to the parcellation table file (for analyses using parcel-based connectivity matrix data)

This will cause the GUI to automatically open the selection windows to the paths/files specified in the default_paths.mat file, which can be convenient if 
e.g., multiple analyses will use the same patient sample, brain mask, etc.

To run the model stacking GUI, type run_stacking_gui into the MATLAB command window while in the main toolkit directory or after adding it to your path.

***Tutorial Notebooks***

Tutorial notebooks can be found in the /tutorials folder. There are tutorials for each type of analysis implemented in the toolkit.

Joseph Griffis, 2024