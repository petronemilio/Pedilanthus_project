import numpy as np
import matplotlib.pyplot as plt
from PIL import Image, ImageOps
#%matplotlib inline
import cv2 as cv
import imutils

#Cargar dos imágenes. Una con la cámara del calabozo y la otra con la cámara del LANCIS
img_1 = cv.imread('../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/EPM6_088.tif')
img_2 = cv.imread('../Data/Pedilanthus/P_tithymaloides/EPM6_S2_1_Inicio_10X.jpg')
#Para ver las dimensiones en pixeles de las imágenes utilizamos shape 
print(img_1.shape,img_2.shape)
#Con imshow vemos las imágenes
cv.imshow("img_1",img_1)
cv.waitKey(0)
cv.imshow("img_2",img_2)
cv.waitKey(0)
#Crear arrreglos de numpy y plotear los histogramas de distribución de los pixeles
img_1_array = np.asarray(img_1)
img_2_array = np.asarray(img_2)
#
plt.hist(img_1_array.ravel())
plt.hist(img_2_array.ravel())
plt.show()
#cv2 guarda por default en formato bgr, entonces es importante transformar a rgb
img_1=cv.cvtColor(img_1, cv.COLOR_BGR2RGB)
img_2=cv.cvtColor(img_2, cv.COLOR_BGR2RGB)

gray_1 = cv.cvtColor(img_1, cv.COLOR_RGB2GRAY)
gray_2 = cv.cvtColor(img_2, cv.COLOR_RGB2GRAY)
plt.hist(img_1.ravel())
plt.hist(gray_1.ravel())
plt.show()

plt.imshow(gray_1, cmap="gray", vmin=0, vmax=255)
plt.show()
plt.imshow(gray_2, cmap="gray", vmin=0, vmax=255)
plt.show()

# global thresholding
ret1,th1_1 = cv.threshold(gray_1,127,255,cv.THRESH_BINARY)
# Otsu's thresholding
ret2,th1_2 = cv.threshold(gray_1,0,255,cv.THRESH_BINARY+cv.THRESH_OTSU)
# Otsu's thresholding after Gaussian filtering
blur_1 = cv.GaussianBlur(gray_1,(5,5),0)
ret3,th1_3 = cv.threshold(blur_1,0,255,cv.THRESH_BINARY+cv.THRESH_OTSU)
# plot all the images and their histograms
images = [gray_1, 0, th1_1,
          gray_1, 0, th1_2,
          blur_1, 0, th1_3]
titles = ['Original Noisy Image','Histogram','Global Thresholding (v=127)',
          'Original Noisy Image','Histogram',"Otsu's Thresholding",
          'Gaussian filtered Image','Histogram',"Otsu's Thresholding"]
for i in range(3):
    plt.subplot(3,3,i*3+1),plt.imshow(images[i*3],'gray')
    plt.title(titles[i*3]), plt.xticks([]), plt.yticks([])
    plt.subplot(3,3,i*3+2),plt.hist(images[i*3].ravel(),256)
    plt.title(titles[i*3+1]), plt.xticks([]), plt.yticks([])
    plt.subplot(3,3,i*3+3),plt.imshow(images[i*3+2],'gray')
    plt.title(titles[i*3+2]), plt.xticks([]), plt.yticks([])



#Pa la otra
# global thresholding
ret1,th2_1 = cv.threshold(gray_2,127,255,cv.THRESH_BINARY)
# Otsu's thresholding
ret2,th2_2 = cv.threshold(gray_2,0,255,cv.THRESH_BINARY+cv.THRESH_OTSU)
# Otsu's thresholding after Gaussian filtering
blur_2 = cv.GaussianBlur(gray_2,(5,5),0)
ret3,th2_3 = cv.threshold(blur_2,0,255,cv.THRESH_BINARY+cv.THRESH_OTSU)
# plot all the images and their histograms
images = [gray_2, 0, th2_1,
          gray_2, 0, th2_2,
          blur_2, 0, th2_3]
titles = ['Original Noisy Image','Histogram','Global Thresholding (v=127)',
          'Original Noisy Image','Histogram',"Otsu's Thresholding",
          'Gaussian filtered Image','Histogram',"Otsu's Thresholding"]
for i in range(3):
    plt.subplot(3,3,i*3+1),plt.imshow(images[i*3],'gray')
    plt.title(titles[i*3]), plt.xticks([]), plt.yticks([])
    plt.subplot(3,3,i*3+2),plt.hist(images[i*3].ravel(),256)
    plt.title(titles[i*3+1]), plt.xticks([]), plt.yticks([])
    plt.subplot(3,3,i*3+3),plt.imshow(images[i*3+2],'gray')
    plt.title(titles[i*3+2]), plt.xticks([]), plt.yticks([])
plt.show()

th2_3 = cv.bitwise_not(th2_3)
plt.imshow(th2_3, vmin=0, vmax=255)
plt.show()

# find contours in the thresholded image
cnts,hierachy=cv.findContours(th2_3,cv.RETR_TREE,cv.CHAIN_APPROX_SIMPLE)

# loop over the contours
for c in cnts:
    # compute the center of the contour
    M = cv.moments(c)
    if M["m00"] != 0:
        cX = int(M["m10"] / M["m00"])
        cY = int(M["m01"] / M["m00"])
    else:
        cX, cY = 0, 0
    # draw the contour and center of the shape on the image
    cv.drawContours(th2_3, [c], -1, (0, 255, 0), 2)
    cv.circle(th2_3, (cX, cY), 7, (255, 255, 255), -1)
    #cv.putText(th2, "center", (cX - 20, cY - 20),cv.FONT_HERSHEY_SIMPLEX, 0.5, (255, 255, 255), 2)
    # display the image
    #cv.imshow("Image", th2)
plt.imshow(th2_3,cmap="gray", vmin=0, vmax=255)
plt.show()

cv.imwrite('centroids.png',th2_3)
#r,g,b = cv.split(img_cv2)
#plt.hist(img_cv2.ravel())

