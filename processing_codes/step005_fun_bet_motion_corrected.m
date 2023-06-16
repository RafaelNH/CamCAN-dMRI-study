function fun_bet_motion_corrected(sub)
cmdbetdwi=['/imaging/local/software/fsl/v5.0.4/x86_64/fsl/bin/bet ', sub, 'DWI_motion ', sub, 'DWI_motion_brain -F -f 0.1 -g 0 -m'];
system(cmdbetdwi);
cmdguzipmsk=['gunzip ', sub, '/DWI_motion_brain_mask.nii.gz'];
system(cmdguzipmsk);
cmdguzipbrain=['gunzip ', sub, '/DWI_motion_brain.nii.gz'];
system(cmdguzipbrain);