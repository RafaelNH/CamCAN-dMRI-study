# -*- coding: utf-8 -*-
"""
Created on Wed Jan 14 13:04:36 2015

@author: rh04

arg1 = path
e.g. Z:\DWI\Subject1\Rep1\

arg2 = name of DWI volume
e.g. DWI_brain
"""

import nibabel as nib
from gibbs_removal import volume_gibbs_removal

from sys import argv

script, arg1, arg2 = argv

# load data
print 'loading ...'
print arg1+arg2+'.nii'
img = nib.load(arg1+arg2+'.nii')
data = img.get_data()

# correct data
print 'Correcting data'
datac, tv = volume_gibbs_removal(data, fn=0, nn=3, fourier=False)

# save data
print 'Saving Data ...'
outpath = arg1+arg2+'_gibbs.nii'
hdr = img.get_header()
affine = hdr.get_base_affine()
im = nib.Nifti1Image(datac, affine)
im.to_filename(outpath)

outpath = arg1+arg2+'_tv.nii'
hdr = img.get_header()
affine = hdr.get_base_affine()
im = nib.Nifti1Image(tv, affine)
im.to_filename(outpath)

print 'Done ...'