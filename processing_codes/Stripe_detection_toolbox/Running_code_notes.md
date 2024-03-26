- The function **fun_white_stripes** takes the paths of the dMRI data and its b-value (for a given subject) as input and computes the apparent 'stripe' index (**aSI**) (Phase 1 described in Appendix A of the paper's supplementary material, cf. Supplementary Fig. S1). You will obtain an aSI value for each DWI volume of that subject.
- After obtaining the apparent 'stripe' index for all subjects, run **full_tripes_index** to compute the full corrected stripe index (**fSI**) (Phase 2 described in Appendix A of the paper's supplementary material, cf. Supplementary Fig. S2). This function takes as input a matrix SI [N_subjects x N_dwi_volumes] containing the aSI values for all subjects and a text file with the b-values of the acquired data (i.e., you only need to specify a b-value file for one subject). It will output the full corrected stripe index [N_subjects x N_dwi_volumes] for all data volumes and all subjects. From visual inspections, this measure has been able to distinguish volumes that were corrupted by motion artifacts from those that were not by setting a threshold value of 0.1. For example, the number of corrupted volumes for a given subject was calculated by counting the number of DWI volumes with an fSI > 0.1.1.