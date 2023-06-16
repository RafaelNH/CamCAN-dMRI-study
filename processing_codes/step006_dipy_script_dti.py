# -*- coding: utf-8 -*-
"""
Created on Wed Jan 14 13:04:36 2015

@author: Rafael Neto Henriques

arg1 = path
e.g. Z:/DWI/Subject1/Rep1/

arg2 = name of DWI in nifti format
e.g. DWI_brain.nii

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
import dipy.reconst.dti as dti
import numpy as np

from sys import argv
from dipy.io import read_bvals_bvecs
from dipy.core.gradients import gradient_table

script, arg1, arg2, arg3, arg4, arg5 = argv

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
print 'fiting DTI ...'
dm = dti.TensorModel(gtab, 'NLLS')
dtifit = dm.fit(data, mask=mask)
FA = dtifit.fa
MD = dtifit.md
RD = dtifit.rd
AD = dtifit.ad
V1 = np.squeeze(dtifit.directions)

# save data
print 'Saving Data ...'
outpath = arg1+'dti/'
if not os.path.isdir(outpath):
    os.mkdir(outpath)
hdr = img.get_header()
affine = hdr.get_base_affine()
imfa = nib.Nifti1Image(FA, affine)
imfa.to_filename(outpath+'FA.nii')
immd = nib.Nifti1Image(MD, affine)
immd.to_filename(outpath+'MD.nii')
imrd = nib.Nifti1Image(RD, affine)
imrd.to_filename(outpath+'RD.nii')
imad = nib.Nifti1Image(AD, affine)
imad.to_filename(outpath+'AD.nii')
imv1 = nib.Nifti1Image(V1, affine)
imv1.to_filename(outpath+'V1.nii')

print 'Done ...'