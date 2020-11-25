#!/usr/bin/python3

import cv2 as cv
import numpy as np
import imutils
import matplotlib.pyplot as plt
import pandas as pd
import filters
import glob
import os
from os import listdir

# get the filenames in directory
epm6_006 = cv.imread("../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/first_series/EPM6_006.tif")

filters06 = filters.getfilters(epm6_006)

#
gray06 = filters.array_gray(filters06)

threshold06 = {}

for k, v in gray06.items():
    threshold06[k] = filters.thresholding(v)
#
threshold_dict = {}
#
threshold_dict.update({'labcol_06':threshold06['Labcolor']['th2'],'blur_06':threshold06['Blur']['th2'],'bilat_06':threshold06['bilateral']['th2']})
#
cv.imwrite('../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/first_series/EPM6_006_labcolor.png',threshold_dict['labcol_06'])
cv.imwrite('../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/first_series/EPM6_006_blur.png',threshold_dict['blur_06'])

centroids = {}
for k, v in threshold_dict.items():
    centroids[k] = filters.out_seg(v)

images_centroids = {}
for k, v in threshold_dict.items():
    images_centroids[k] = filters.out_seg_fig(v)

cv.imwrite('../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/first_series/EPM6_006_labcolor_centr.png',images_centroids['labcol_06'])
cv.imwrite('../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/first_series/EPM6_006_blur_centr.png',images_centroids['blur_06'])
cv.imwrite('../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/first_series/EPM6_006_bilat_centr.png',images_centroids['bilat_06'])

##################3
# get the filenames in directory
epm6_007 = cv.imread("../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/first_series/EPM6_007.tif")

filters07 = filters.getfilters(epm6_007)

#
gray07 = filters.array_gray(filters07)

threshold07 = {}

for k, v in gray07.items():
    threshold07[k] = filters.thresholding(v)
#
threshold_dict = {}
#
threshold_dict.update({'labcol_07':threshold07['Labcolor']['th2'],'blur_07':threshold07['Blur']['th2'],'bilat_07':threshold07['bilateral']['th2']})

cv.imwrite('../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/first_series/EPM6_007_labcolor.png',threshold_dict['labcol_07'])
cv.imwrite('../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/first_series/EPM6_007_blur.png',threshold_dict['blur_07'])
#####
centroids = {}
for k, v in threshold_dict.items():
    centroids[k] = filters.out_seg(v)

images_centroids = {}
for k, v in threshold_dict.items():
    images_centroids[k] = filters.out_seg_fig(v)

cv.imwrite('../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/first_series/EPM6_007_labcolor_centr.png',images_centroids['labcol_07'])
cv.imwrite('../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/first_series/EPM6_007_blur_centr.png',images_centroids['blur_07'])
cv.imwrite('../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/first_series/EPM6_007_bilat_centr.png',images_centroids['bilat_07'])

#####################################
# get the filenames in directory
epm6_008 = cv.imread("../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/first_series/EPM6_008.tif")

filters08 = filters.getfilters(epm6_008)

#
gray08 = filters.array_gray(filters08)

threshold08 = {}

for k, v in gray08.items():
    threshold08[k] = filters.thresholding(v)
#
threshold_dict = {}
#
threshold_dict.update({'labcol_08':threshold08['Labcolor']['th2'],'blur_08':threshold08['Blur']['th2'],'bilat_08':threshold08['bilateral']['th2']})

cv.imwrite('../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/first_series/EPM6_008_labcolor.png',threshold_dict['labcol_08'])
cv.imwrite('../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/first_series/EPM6_008_blur.png',threshold_dict['blur_08'])

centroids = {}
for k, v in threshold_dict.items():
    centroids[k] = filters.out_seg(v)

images_centroids = {}
for k, v in threshold_dict.items():
    images_centroids[k] = filters.out_seg_fig(v)

cv.imwrite('../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/first_series/EPM6_008_labcolor_centr.png',images_centroids['labcol_08'])
cv.imwrite('../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/first_series/EPM6_008_blur_centr.png',images_centroids['blur_08'])
cv.imwrite('../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/first_series/EPM6_008_bilat_centr.png',images_centroids['bilat_08'])

##########################
# get the filenames in directory
epm6_009 = cv.imread("../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/first_series/EPM6_009.tif")

filters09 = filters.getfilters(epm6_009)

#
gray09 = filters.array_gray(filters09)

threshold09 = {}

for k, v in gray09.items():
    threshold09[k] = filters.thresholding(v)
#
threshold_dict = {}
#
threshold_dict.update({'labcol_09':threshold09['Labcolor']['th2'],'blur_09':threshold09['Blur']['th2'],'bilat_09':threshold09['bilateral']['th2']})

cv.imwrite('../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/first_series/EPM6_009_labcolor.png',threshold_dict['labcol_09'])
cv.imwrite('../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/first_series/EPM6_009_blur.png',threshold_dict['blur_09'])

centroids = {}
for k, v in threshold_dict.items():
    centroids[k] = filters.out_seg(v)

images_centroids = {}
for k, v in threshold_dict.items():
    images_centroids[k] = filters.out_seg_fig(v)

cv.imwrite('../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/first_series/EPM6_009_labcolor_centr.png',images_centroids['labcol_09'])
cv.imwrite('../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/first_series/EPM6_009_blur_centr.png',images_centroids['blur_09'])
cv.imwrite('../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/first_series/EPM6_009_bilat_centr.png',images_centroids['bilat_09'])
#################################
# get the filenames in directory
epm6_010 = cv.imread("../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/first_series/EPM6_010.tif")

filters10 = filters.getfilters(epm6_010)

#
gray10 = filters.array_gray(filters10)

threshold10 = {}

for k, v in gray10.items():
    threshold10[k] = filters.thresholding(v)
#
threshold_dict = {}
#
threshold_dict.update({'labcol_10':threshold10['Labcolor']['th2'],'blur_10':threshold10['Blur']['th2'],'bilat_10':threshold10['bilateral']['th2']})
###
cv.imwrite('../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/first_series/EPM6_010_labcolor.png',threshold_dict['labcol_10'])
cv.imwrite('../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/first_series/EPM6_010_blur.png',threshold_dict['blur_10'])
#
centroids = {}
for k, v in threshold_dict.items():
    centroids[k] = filters.out_seg(v)

images_centroids = {}
for k, v in threshold_dict.items():
    images_centroids[k] = filters.out_seg_fig(v)

cv.imwrite('../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/first_series/EPM6_010_labcolor_centr.png',images_centroids['labcol_10'])
cv.imwrite('../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/first_series/EPM6_010_blur_centr.png',images_centroids['blur_10'])
cv.imwrite('../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/first_series/EPM6_010_bilat_centr.png',images_centroids['bilat_10'])

############# get the filenames in directory
epm6_011 = cv.imread("../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/first_series/EPM6_011.tif")

filters11 = filters.getfilters(epm6_011)

#
gray11 = filters.array_gray(filters11)

threshold11 = {}

for k, v in gray11.items():
    threshold11[k] = filters.thresholding(v)
#
threshold_dict = {}
#
threshold_dict.update({'labcol_11':threshold11['Labcolor']['th2'],'blur_11':threshold11['Blur']['th2'],'bilat_11':threshold11['bilateral']['th2']})

cv.imwrite('../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/first_series/EPM6_011_labcolor.png',threshold_dict['labcol_11'])
cv.imwrite('../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/first_series/EPM6_011_blur.png',threshold_dict['blur_11'])

centroids = {}
for k, v in threshold_dict.items():
    centroids[k] = filters.out_seg(v)

images_centroids = {}
for k, v in threshold_dict.items():
    images_centroids[k] = filters.out_seg_fig(v)

cv.imwrite('../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/first_series/EPM6_011_labcolor_centr.png',images_centroids['labcol_11'])
cv.imwrite('../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/first_series/EPM6_011_blur_centr.png',images_centroids['blur_11'])
cv.imwrite('../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/first_series/EPM6_011_bilat_centr.png',images_centroids['bilat_11'])


thin06 = cv.ximgproc.thinning(threshold06['Labcolor']['th2'])
cv.imwrite('../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/first_series/EPM6_006_thin.png', thin06)

thin07 = cv.ximgproc.thinning(threshold07['Labcolor']['th2'])
cv.imwrite('../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/first_series/EPM6_007_thin.png', thin07)

thin08 = cv.ximgproc.thinning(threshold08['Labcolor']['th2'])
cv.imwrite('../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/first_series/EPM6_008_thin.png', thin08)

thin09 = cv.ximgproc.thinning(threshold09['Labcolor']['th2'])
cv.imwrite('../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/first_series/EPM6_009_thin.png', thin09)

thin10 = cv.ximgproc.thinning(threshold10['Labcolor']['th2'])
cv.imwrite('../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/first_series/EPM6_010_thin.png', thin10)

thin11 = cv.ximgproc.thinning(threshold11['Labcolor']['th2'])
cv.imwrite('../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/first_series/EPM6_011_thin.png', thin11)




