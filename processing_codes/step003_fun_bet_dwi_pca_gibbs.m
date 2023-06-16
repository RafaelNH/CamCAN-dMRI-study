function fun_bet_dwi_pca_gibbs(sub)
cmdbetdwi=['/imaging/local/software/fsl/v5.0.4/x86_64/fsl/bin/bet ', sub, 'DWI_pca_gibbs ', sub, 'DWI_pca_gibbs_brain -F -f 0.1 -g 0 -m'];
system(cmdbetdwi);
cmdguzipmsk=['gunzip ', sub, '/DWI_pca_gibbs_brain_mask.nii.gz'];
system(cmdguzipmsk);
cmdguzipbrain=['gunzip ', sub, '/DWI_pca_gibbs_brain.nii.gz'];
system(cmdguzipbrain);