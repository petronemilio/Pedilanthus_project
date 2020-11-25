#!/usr/bin/python3

import cv2 as cv
import numpy as np
import imutils
import matplotlib.pyplot as plt
import pandas as pd
import filters

######
epm6_011 = cv.imread('../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/EPM6_011.tif')
epm6_012 = cv.imread('../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/EPM6_012.tif')
epm6_013 = cv.imread('../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/EPM6_013.tif')
epm6_014 = cv.imread('../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/EPM6_014.tif')
epm6_015 = cv.imread('../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/EPM6_015.tif')

filters11 = filters.getfilters(epm6_011)
filters12 = filters.getfilters(epm6_012)
filters13 = filters.getfilters(epm6_013)
filters14 = filters.getfilters(epm6_014)
filters15 = filters.getfilters(epm6_015)

#
gray11 = filters.array_gray(filters11)
gray12 = filters.array_gray(filters12)
gray13 = filters.array_gray(filters13)
gray14 = filters.array_gray(filters14)
gray15 = filters.array_gray(filters15)


threshold11 = {}
threshold12 = {}
threshold13 = {}
threshold14 = {}
threshold15 = {}

for k, v in gray11.items():
    threshold11[k] = filters.thresholding(v)

cv.imwrite('../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/prueba/EPM6_011_labcolorth1.png',threshold11['Labcolor']['th1'])
cv.imwrite('../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/prueba/EPM6_011_labcolorth2.png',threshold11['Labcolor']['th2'])
cv.imwrite('../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/prueba/EPM6_011_distanceth1.png',threshold11['Distance']['th1'])
cv.imwrite('../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/prueba/EPM6_011_distanceth2.png',threshold11['Distance']['th2'])
cv.imwrite('../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/prueba/EPM6_011_blurth1.png',threshold11['Blur']['th1'])
cv.imwrite('../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/prueba/EPM6_011_blurth2.png',threshold11['Blur']['th2'])
cv.imwrite('../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/prueba/EPM6_011_medianth1.png',threshold11['median']['th1'])
cv.imwrite('../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/prueba/EPM6_011_medianth2.png',threshold11['median']['th2'])
cv.imwrite('../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/prueba/EPM6_011_bilateralth1.png',threshold11['bilateral']['th1'])
cv.imwrite('../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/prueba/EPM6_011_bilateralth2.png',threshold11['bilateral']['th2'])
cv.imwrite('../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/prueba/EPM6_011_shiftedth1.png',threshold11['shifted']['th1'])
cv.imwrite('../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/prueba/EPM6_011_shiftedth2.png',threshold11['shifted']['th2'])

for k, v in gray12.items():
    threshold12[k] = filters.thresholding(v)

cv.imwrite('../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/prueba/EPM6_012_labcolorth1.png',threshold12['Labcolor']['th1'])
cv.imwrite('../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/prueba/EPM6_012_labcolorth2.png',threshold12['Labcolor']['th2'])
cv.imwrite('../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/prueba/EPM6_012_distanceth1.png',threshold12['Distance']['th1'])
cv.imwrite('../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/prueba/EPM6_012_distanceth2.png',threshold12['Distance']['th2'])
cv.imwrite('../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/prueba/EPM6_012_blurth1.png',threshold12['Blur']['th1'])
cv.imwrite('../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/prueba/EPM6_012_blurth2.png',threshold12['Blur']['th2'])
cv.imwrite('../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/prueba/EPM6_012_medianth1.png',threshold12['median']['th1'])
cv.imwrite('../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/prueba/EPM6_012_medianth2.png',threshold12['median']['th2'])
cv.imwrite('../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/prueba/EPM6_012_bilateralth1.png',threshold12['bilateral']['th1'])
cv.imwrite('../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/prueba/EPM6_012_bilateralth2.png',threshold12['bilateral']['th2'])
cv.imwrite('../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/prueba/EPM6_012_shiftedth1.png',threshold12['shifted']['th1'])
cv.imwrite('../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/prueba/EPM6_012_shiftedth2.png',threshold12['shifted']['th2'])

for k, v in gray13.items():
    threshold13[k] = filters.thresholding(v)

cv.imwrite('../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/prueba/EPM6_013_labcolorth1.png',threshold13['Labcolor']['th1'])
cv.imwrite('../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/prueba/EPM6_013_labcolorth2.png',threshold13['Labcolor']['th2'])
cv.imwrite('../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/prueba/EPM6_013_distanceth1.png',threshold13['Distance']['th1'])
cv.imwrite('../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/prueba/EPM6_013_distanceth2.png',threshold13['Distance']['th2'])
cv.imwrite('../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/prueba/EPM6_013_blurth1.png',threshold13['Blur']['th1'])
cv.imwrite('../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/prueba/EPM6_013_blurth2.png',threshold13['Blur']['th2'])
cv.imwrite('../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/prueba/EPM6_013_medianth1.png',threshold13['median']['th1'])
cv.imwrite('../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/prueba/EPM6_013_medianth2.png',threshold13['median']['th2'])
cv.imwrite('../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/prueba/EPM6_013_bilateralth1.png',threshold13['bilateral']['th1'])
cv.imwrite('../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/prueba/EPM6_013_bilateralth2.png',threshold13['bilateral']['th2'])
cv.imwrite('../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/prueba/EPM6_013_shiftedth1.png',threshold13['shifted']['th1'])
cv.imwrite('../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/prueba/EPM6_013_shiftedth2.png',threshold13['shifted']['th2'])

for k, v in gray14.items():
    threshold14[k] = filters.thresholding(v)

cv.imwrite('../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/prueba/EPM6_014_labcolorth1.png',threshold14['Labcolor']['th1'])
cv.imwrite('../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/prueba/EPM6_014_labcolorth2.png',threshold14['Labcolor']['th2'])
cv.imwrite('../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/prueba/EPM6_014_distanceth1.png',threshold14['Distance']['th1'])
cv.imwrite('../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/prueba/EPM6_014_distanceth2.png',threshold14['Distance']['th2'])
cv.imwrite('../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/prueba/EPM6_014_blurth1.png',threshold14['Blur']['th1'])
cv.imwrite('../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/prueba/EPM6_014_blurth2.png',threshold14['Blur']['th2'])
cv.imwrite('../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/prueba/EPM6_014_medianth1.png',threshold14['median']['th1'])
cv.imwrite('../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/prueba/EPM6_014_medianth2.png',threshold14['median']['th2'])
cv.imwrite('../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/prueba/EPM6_014_bilateralth1.png',threshold14['bilateral']['th1'])
cv.imwrite('../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/prueba/EPM6_014_bilateralth2.png',threshold14['bilateral']['th2'])
cv.imwrite('../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/prueba/EPM6_014_shiftedth1.png',threshold14['shifted']['th1'])
cv.imwrite('../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/prueba/EPM6_014_shiftedth2.png',threshold14['shifted']['th2'])

for k, v in gray15.items():
    threshold15[k] = filters.thresholding(v)

cv.imwrite('../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/prueba/EPM6_015_labcolorth1.png',threshold15['Labcolor']['th1'])
cv.imwrite('../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/prueba/EPM6_015_labcolorth2.png',threshold15['Labcolor']['th2'])
cv.imwrite('../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/prueba/EPM6_015_distanceth1.png',threshold15['Distance']['th1'])
cv.imwrite('../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/prueba/EPM6_015_distanceth2.png',threshold15['Distance']['th2'])
cv.imwrite('../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/prueba/EPM6_015_blurth1.png',threshold15['Blur']['th1'])
cv.imwrite('../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/prueba/EPM6_015_blurth2.png',threshold15['Blur']['th2'])
cv.imwrite('../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/prueba/EPM6_015_medianth1.png',threshold15['median']['th1'])
cv.imwrite('../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/prueba/EPM6_015_medianth2.png',threshold15['median']['th2'])
cv.imwrite('../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/prueba/EPM6_015_bilateralth1.png',threshold15['bilateral']['th1'])
cv.imwrite('../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/prueba/EPM6_015_bilateralth2.png',threshold15['bilateral']['th2'])
cv.imwrite('../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/prueba/EPM6_015_shiftedth1.png',threshold15['shifted']['th1'])
cv.imwrite('../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/prueba/EPM6_015_shiftedth2.png',threshold15['shifted']['th2'])



threshold_dict = {}

threshold_dict.update({'labcol_11':threshold11['Labcolor']['th2'], 'dst_11':threshold11['Distance']['th2'],
                    'blur_11':threshold11['Blur']['th2'],'median_11':threshold11['median']['th2'],
                    'bilat_11':threshold11['bilateral']['th2'], 'shift_11':threshold11['shifted']['th2'],
                    'labcol_12':threshold12['Labcolor']['th2'], 'dst_12':threshold12['Distance']['th2'],
                    'blur_12':threshold12['Blur']['th2'],'median_12':threshold12['median']['th2'],
                    'bilat_12':threshold12['bilateral']['th2'], 'shift_12':threshold12['shifted']['th2'],
                    'labcol_13':threshold13['Labcolor']['th2'], 'dst_13':threshold13['Distance']['th2'],
                    'blur_13':threshold13['Blur']['th2'],'median_13':threshold13['median']['th2'],
                    'bilat_13':threshold13['bilateral']['th2'], 'shift_13':threshold13['shifted']['th2'],
                    'labcol_14':threshold14['Labcolor']['th2'], 'dst_14':threshold14['Distance']['th2'],
                    'blur_14':threshold14['Blur']['th2'],'median_14':threshold14['median']['th2'],
                    'bilat_14':threshold14['bilateral']['th2'], 'shift_14':threshold14['shifted']['th2'],
                    'labcol_15':threshold15['Labcolor']['th2'], 'dst_15':threshold15['Distance']['th2'],
                    'blur_15':threshold15['Blur']['th2'],'median_15':threshold15['median']['th2'],
                    'bilat_15':threshold15['bilateral']['th2'], 'shift_15':threshold15['shifted']['th2']})
##Calcular centroides y otros parámetros 
centroids = {}
for k, v in threshold_dict.items():
    centroids[k] = filters.out_seg(v)


#centroids_df = pd.DataFrame(centroids)

filter_ids = centroids.keys()
appended_data = []

for keys in filter_ids:
    temp_df = pd.DataFrame.from_dict(centroids[keys])
    temp_df.insert(0,"Segmentation", np.repeat(keys, np.shape(temp_df)[0]), True)
    appended_data.append(temp_df)
    
appended_data = pd.concat(appended_data) 

appended_data.to_csv('../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/prueba/parameters_segmentation.csv', index=False)


#centroids_df.to_csv('../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/prueba/centroids.csv', index=False)


