# -*- coding: utf-8 -*-
"""
Created on Mon Oct 11 10:20:52 2021

@author: nilou
"""

import sys
sys.path.append('')

import io
import os
import cv2 as cv
import skimage
import glob
import csv
import math
import numpy as np
import scipy as sp
import pandas as pd
import matplotlib 
import napari
from numpy.core.umath_tests import inner1d
from matplotlib import pyplot as plt
from scipy import ndimage, misc
from skimage import io, color, measure, segmentation
from skimage.measure import label, regionprops, regionprops_table
import czifile as zis
from aicsimageio import AICSImage
from aicsimageio.writers import OmeTiffWriter
from scipy import ndimage
from skimage.morphology import erosion, dilation, opening, closing, white_tophat, disk, black_tophat, skeletonize, convex_hull_image
import cv2 as cv
import matplotlib.pyplot as plt
from pathlib import Path
from skimage.filters import threshold_otsu, threshold_local


# 1: reading and loading the image from directory ######################################

PATH_in = ''
PATH_out = ''     #### specify the path to save the images

name = Path(PATH_in).stem
pixel_to_um = 1                                             #### specify scale of the image
image = AICSImage(PATH_in)
# img = czifile.imread(PATH_in)
ch = image.get_image_data("YX", S=1, T=0 , Z=0, C=0)  # returns 3D ZYX numpy array 6D STCZYX (.data)
ch1 = image.get_image_data("YX", S=1, T=0 , Z=0, C=2) 
Channels = image.channel_names

print(np.max(ch), np.max(ch1), Channels)
ch = ch/np.max(ch)
ch1 = ch1/np.max(ch1)



#%%##### REFMASK #################################################################################################
# REFMASK FUNCTION##############################################################################################
def ref_mask(channel_index):   
    channel = image.get_image_data("YX", C=channel_index)
    
    ref_mask = threshold_otsu(channel)
    ref_mask = channel > np.average(ref_mask)*0.1
    # # for i in range(1):
    # # ref_mask = ndimage.binary_opening(ref_mask, iterations=1)
    # # ref_mask = ndimage.binary_closing(ref_mask, iterations=1)
    ref_mask = ndimage.binary_erosion(ref_mask, iterations=1)
    # # ref_mask = ndimage.binary_fill_holes(ref_mask)
    
    return(ref_mask)
##################################################################################################

#%%

# ref_mask_int = ref_mask.astype(np.int) 
ch_dox = ch * ref_mask(0)
# ch_dox = (ch_dox - np.min(ch_dox)) / (np.max(ch_dox) - np.min(ch_dox))
ch_particles = ch1 * ref_mask(2)
# ch_particles = (ch_particles - np.min(ch_particles)) / (np.max(ch_particles) - np.min(ch_particles))


fig, ax = plt.subplots(1,2, subplot_kw=dict(aspect="equal"), figsize=[25, 20])
ax[0].imshow(ch_dox>0)
ax[0].set_title('Doxorubicin '+ name, size=30)
ax[0].axis('off')

ax[1].imshow(ch_particles>0)
ax[1].set_title('Particles ' + name, size=30)
ax[1].axis('off')

plt.savefig (PATH_out + name + '_corr.png', bbox_inches='tight')
plt.tight_layout()
plt.show()

    
#%%    # ### making Fluorograms vectors of intensities #################################################
# STATISTICS #####################################################################################
### making a 2d array to save vectors of each channel Intensities and mask #######################

corr = pd.DataFrame()

dox_i = ch_dox.ravel()
particle_i = ch_particles.ravel()
REF = ref_mask(0).ravel()

dox_i_pos = dox_i[REF>0]
particle_i_pos = particle_i[REF>0]

corr ['Filename'] = [name for x in range(len(dox_i_pos))]
corr ['dox'] = dox_i_pos [:]
corr ['Particles'] = particle_i_pos [:]


# split k-overlap coefficients

import scipy.stats as stats

k_dox = np.dot(dox_i_pos, particle_i_pos)/np.dot(dox_i_pos, dox_i_pos)
corr["k_dox"] = k_dox
k_particles = np.dot(dox_i_pos, particle_i_pos)/np.dot(particle_i_pos, particle_i_pos)
corr["k_particles"] = k_particles
Manders = math.sqrt (k_dox * k_particles)
corr["Manders"] = Manders
r, p = stats.pearsonr(dox_i_pos, particle_i_pos)
corr["Pearson"] = r
corr["P"] = p
tau, pk = stats.kendalltau(dox_i_pos, particle_i_pos)
corr["Kendall"] = tau
corr["Pk"] = pk
rho, ps = stats.spearmanr(dox_i_pos, particle_i_pos)
corr["Spearman"] = rho
corr["Ps"] = ps

#### Writing STATISTICS File ######################################################################

corr.to_csv(PATH_out + name + '_corr.csv', index=False)

#%%#######################################################################################################
# Calculate the point density

# from scipy.stats import gaussian_kde
# xy = np.vstack([dox_i,particle_i])
# z = gaussian_kde(xy)(xy)

import seaborn as sns
j = sns.jointplot(x= dox_i_pos, y= particle_i_pos, xlim=[0, 1], ylim=[0, 1],
                  kind='reg', scatter_kws={"s": 10, "alpha": 0.15})

j.ax_joint.legend(["Pearson ={:.2e}, p = {:.2e}".format(r,p),
                   "Manders ={:.2e}".format(Manders)])
j.ax_joint.set_xlabel('Doxorubicin', fontsize=14)
j.ax_joint.set_ylabel('Nanoparticles', fontsize=14)
# j.ax_joint.set_title(name, fontsize=14, loc='center')
plt.savefig (PATH_out + name + '_jointplot.png', bbox_inches='tight')
plt.show()

    
##################################################################################################











