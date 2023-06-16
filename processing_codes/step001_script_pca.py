# -*- coding: utf-8 -*-
"""
Created on 06/03/2017

@author: rh04

arg1 = path
e.g. Z:\DWI\Subject1\Rep1\

arg2 = name of DWI volume
e.g. DWI_brain
"""

import nibabel as nib
import numpy as np
from pca_utils import pca_denoising

from sys import argv

script, arg1, arg2 = argv

# load data
print 'loading ...'
print arg1+arg2+'.nii'
img = nib.load(arg1+arg2+'.nii')
data = img.get_data()

# correct data
print 'Correcting data'
datac, std, ncomps = pca_denoising(data, ps=2, overcomplete=True)

affine = np.identity(4)

# save data
print 'Saving Data ...'
outpath = arg1+arg2+'_pca.nii'
im = nib.Nifti1Image(datac, affine)
im.to_filename(outpath)

outpath = arg1+arg2+'_std.nii'
im = nib.Nifti1Image(std, affine)
im.to_filename(outpath)

outpath = arg1+arg2+'_ncomps.nii'
im = nib.Nifti1Image(ncomps, affine)
im.to_filename(outpath)

print 'Done ...'