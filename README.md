# CamCAN-dMRI-study
 
This repository contains the code used for the following research paper:

Henriques, RN, Henson, R, Cam-CAN, Correia, MM. Unique information from common diffusion MRI models about white-matter differences across the human adult lifespan. 2023 (currently under revision)

# Content

- Folder "processing_codes": contains the general code used to process DTI/DKI/NODDI. Note, only the essential steps to allow the replication of our data processing were posted. Cam-CAN raw data is available at https://camcan-archive.mrc-cbu.cam.ac.uk/dataaccess/. If you want to run the processing codes after having access to the Cam-CAN data you will have to insert our own data paths in the scripts. Some code packages are also missing (e.g. NODDI toolbox, NIFTI toolbox, DIPY python codes) which can be found on-line.

- Folder "data": Containing a matlab file with the final output from the last processing step in "processing_codes". This allows running the analysing scripts without the need of reruning the code in "processing_codes".

- Folder "analysis_codes": contains all the code to generate all the figures of the paper. You will be able to run the main script "main_dMRI_stats_analysis.m" from the data saved in folder "data", after downloading the following from MathWorks file exchange: fdr_bh.m (https://www.mathworks.com/matlabcentral/fileexchange/27418-fdr_bh), redblue.m (https://it.mathworks.com/matlabcentral/fileexchange/25536-red-blue-colormap), turbo.m (https://it.mathworks.com/matlabcentral/fileexchange/74662-turbo).

Note: If you are only interested to have the final diffusion metric values averaged across WM voxels and averaged separetely for individual WM ROIs, these are available in the CSV files “Global_Metrics.csv” (averaged across WM voxels) and “ROI_Metrics.csv” (separately for each ROI).

