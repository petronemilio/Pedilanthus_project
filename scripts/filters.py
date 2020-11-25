#!/usr/bin/python3

import cv2 as cv
import numpy as np

def labcolor(image):
    """ Labcolor is a function to transform the image with the lab method"""
    rgb_img = cv.cvtColor(image, cv.COLOR_BGR2RGB)
    lab = cv.cvtColor(rgb_img, cv.COLOR_RGB2LAB)
    l,a,b = cv.split(lab)
    clahe = cv.createCLAHE(clipLimit = 3.0, tileGridSize = (8,8))
    cl = clahe.apply(l)
    limg = cv.merge((cl, a, b))
    lab_rgb = cv.cvtColor(limg, cv.COLOR_LAB2RGB)
    return lab_rgb 

def getfilters(image):
    """ This function retrieves an image in different filters.
        Avoids writing repetive tasks. """
    images = {}
    rgb_img = cv.cvtColor(image, cv.COLOR_BGR2RGB)
    #Add labcolor to images
    lab_rgb = labcolor(image)
    images["Labcolor"] = lab_rgb
    #Distance filter
    kernel = np.ones((7,7), np.float32)/49
    dst = cv.filter2D(rgb_img, -1, kernel)
    images["Distance"] = dst
    #Blur filter
    blur = cv.GaussianBlur(rgb_img, (5,5), 0)
    images["Blur"] = blur
    #Median filter
    median = cv.medianBlur(rgb_img, 5)
    images["median"] = median   
    #Bilateral filter
    bil = cv.bilateralFilter(rgb_img, 9,200,200)
    images["bilateral"] = bil
    #Shifted
    shifted = cv.pyrMeanShiftFiltering(rgb_img, 21, 51)
    images["shifted"] = shifted
    return images

def array_gray(dictionary):
    """ Recieves a dictionary with an array of images an gives back
        an array of gray images"""
    gray={}
    for k, v in dictionary.items():
        x = cv.cvtColor(v,cv.COLOR_BGR2GRAY)
        gray[k] = cv.equalizeHist(x)
    return gray

def thresholding(image):
    """The function recieves an image an performs to types of thresholding:
     global and Otsu's""" 
    threshold = {}
    # global thresholing
    ret1,th1 = cv.threshold(image,127,255,cv.THRESH_BINARY) #dst filter
    # Otsu's thresholding
    ret2,th2 = cv.threshold(image,0,255,cv.THRESH_BINARY+cv.THRESH_OTSU) #dst filter
    th1_inv = cv.bitwise_not(th1)
    th2_inv = cv.bitwise_not(th2)
    threshold["th1"] = th1_inv
    threshold["th2"] = th2_inv
    return threshold


def out_seg(image):
    # find contours in the thresholded image
    cnts,hierarchy = cv.findContours(image, cv.RETR_TREE,cv.CHAIN_APPROX_NONE)#SIMPLE#RETR_EXTERNAL
    param_img={"centroids":[],"areas":[], "perimeters":[]}
    #make array of centroids
    #Array of areas
    #Array for perimeters
    # loop over the contours
    for c in cnts:
    # compute the center of the contour
        M = cv.moments(c) #que es moments?¡¿
        if M["m00"] != 0:
            cX = int(M["m10"] / M["m00"])
            cY = int(M["m01"] / M["m00"])
        else:
            cX, cY = 0, 0
    # draw the contour and center of the shape on the image
        #cv.drawContours(image, [c], -1, (0, 255, 0), 2)
        #cv.circle(th3_2, (cX, cY), 7, (255, 255, 255), -1)
        param_img["centroids"].append([cX, cY])
        param_img["areas"].append(cv.contourArea(c))
        param_img["perimeters"].append(cv.arcLength(c,True))
    return param_img

def out_seg_fig(image):
    image_copy = image.copy()
    # find contours in the thresholded image
    cnts,hierarchy = cv.findContours(image_copy, cv.RETR_TREE,cv.CHAIN_APPROX_NONE)#SIMPLE#RETR_EXTERNAL
    #make array of centroids
    centroids = []
    #Array of areas
    #Array for perimeters
    # loop over the contours
    for c in cnts:
    # compute the center of the contour
        M = cv.moments(c) #que es moments?¡¿
        if M["m00"] != 0:
            cX = int(M["m10"] / M["m00"])
            cY = int(M["m01"] / M["m00"])
        else:
            cX, cY = 0, 0
        # draw the contour and center of the shape on the image
        cv.drawContours(image_copy, [c], -1, (0, 255, 0), 2)
        #cv.circle(th3_2, (cX, cY), 7, (255, 255, 255), -1)
        centroids.append([cX, cY])
    return image_copy






 
