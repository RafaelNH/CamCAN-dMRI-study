# -*- coding: utf-8 -*-
"""
Created on Wed July 14 13:04:36 2016

@author: rh04

arg1 = path
e.g. Z:/DWI/Subject1/Rep1/

arg2 = name of DWI in nifti format
e.g. DWI.nii

arg3 = name of DWI mask in nifti format
e.g. DWI_brain_mask.nii

arg4 = name of bvals in bval format
e.g. bvals.bval

arg5 = name of bvecs in bvec format
e.g. bvecs.bvec

make sure you have the python path in .cshrc. For this add the following lines
to .cshrc

setenv PYTHONPATH /home/'your_user_name'/
setenv PYTHONPATH /imaging/local/software/python_packages/dipy/0.9.2/:${PYTHONPATH}
setenv PYTHONPATH /imaging/local/software/python_packages/nibabel/1.3.0/:${PYTHONPATH}
"""

import os
import nibabel as nib
import numpy as np
from dki_alternative import avs_dki_df

from sys import argv
from dipy.io import read_bvals_bvecs
from dipy.core.gradients import gradient_table

script, arg1, arg2, arg3, arg4, arg5 = argv

#arg1 = '/imaging/rh04/camcan_dki/DWI05/CBU110220_MR10033_CC610050/'
#arg2 = 'DWI_motion_brain.nii'
#arg3 = 'DWI_motion_brain_mask.nii'
#arg4 = 'bvals.bval'
#arg5 = 'bvecs_orient.bvec'

# load data
print 'loading data ...'
img = nib.load(arg1+arg2)
data = img.get_data()

print 'loading mask ...'
imask = nib.load(arg1+arg3)
mask = imask.get_data()

# load bvals and bvecs and convert them to a gradient table
bvals, bvecs = read_bvals_bvecs(arg1+arg4, arg1+arg5)
gtab = gradient_table(bvals, bvecs)

# fit DTI
print 'fiting DKI ...'
params = avs_dki_df(gtab, data, mask=mask)

print 'extract DKI derived measures ...'
MD = params[..., 0]
MK = params[..., 1]
S0 = params[..., 2]

# save data
print 'Saving Data ...'
outpath = arg1+'dkidf/'
if not os.path.isdir(outpath):
    os.mkdir(outpath)
hdr = img.get_header()
affine = hdr.get_base_affine()
immd = nib.Nifti1Image(MD, affine)
immd.to_filename(outpath+'MD.nii')
immk = nib.Nifti1Image(MK, affine)
immk.to_filename(outpath+'MK.nii')
imrd = nib.Nifti1Image(S0, affine)
imrd.to_filename(outpath+'S0.nii')

print 'Done ...'