#Calculate the centroid of various pixel regions from a CCD image
#06/24/2018
#Dahlia Dry
import numpy as np
import math
from astropy.io import fits
import matplotlib.pyplot as plt

class Centroid(object):
    """
    filepath(str) : filepath of image
    dims(list of length 2) : dimensions of subfield in image you want to centroid
    center(list of length 2): coordinates of center of subfield you want to centroid
    """
    def __init__(self, filepath, dims =[], center = []):
        self.filepath = filepath
        self.img_array = fits.getdata(self.filepath)
        self.dims = dims
        self.center = center

    def display(self):
        plt.imshow(self.img_array,vmin = self.img_array.mean(),
                   vmax = 2*self.img_array.mean())
        plt.gray()
        plt.show()

    def get_sums(self, img=True, array=[]):
        x_sums = []
        y_sums = []
        placeholder = 0
        if img and len(self.dims) == 0 and len(self.center) == 0 :
            for i in self.img_array:
                for j in i:
                    placeholder = placeholder + j
                y_sums.append(placeholder)
                placeholder = 0
            placeholder = 0
            for i in range(len(self.img_array[0])):
                for j in self.img_array[:,i]:
                    placeholder = placeholder + j
                x_sums.append(placeholder)
                placeholder = 0

        elif img and len(self.center) != 0 and len(self.dims) == 0:
            pix = self.img_array[self.center[1]-1:self.center[1]+2,
                                 self.center[0]-1:self.center[0]+2]
            pix = np.asarray(pix)
            for i in pix:
                for j in i:
                    placeholder = placeholder + j
                y_sums.append(placeholder)
                placeholder = 0
            placeholder = 0
            for i in range(len(pix[0])):
                for j in pix[:,i]:
                    placeholder = placeholder + j
                x_sums.append(placeholder)
                placeholder = 0

        elif img and len(self.center) != 0 and len(self.dims) != 0:
            pix = self.img_array[self.center[1]-math.floor(self.dims[1]/2):
                                 self.center[1]+math.floor(self.dims[1]/2)+1,
                                 self.center[0]-math.floor(self.dims[0]/2):
                                 self.center[0]+math.floor(self.dims[0]/2)+1]
            pix = np.asarray(pix)
            for i in pix:
                for j in i:
                    placeholder = placeholder + j
                y_sums.append(placeholder)
                placeholder = 0
            placeholder = 0
            for i in range(len(pix[0])):
                for j in pix[:,i]:
                    placeholder = placeholder + j
                x_sums.append(placeholder)
                placeholder = 0

        elif len(array) != 0:
            for i in array:
                for j in i:
                    placeholder = placeholder + j
                y_sums.append(placeholder)
                placeholder = 0
            placeholder = 0
            for i in range(len(array[0])):
                for j in array[:,i]:
                    placeholder = placeholder + j
                x_sums.append(placeholder)
                placeholder = 0
        else:
            print('no condition satisfied')
        return x_sums,y_sums

    def get_centroids(self, img=True, array = []):
        if img and len(self.center) != 0 and len(self.dims) == 0:
            x_sums, y_sums = self.get_sums()
            x_indices = [self.center[0]-1, self.center[0], self.center[0]+1]
            y_indices = [self.center[1]-1, self.center[1], self.center[1]+1]

        elif img and len(self.center) != 0 and len(self.dims) != 0:
            x_sums, y_sums = self.get_sums()
            x_indices = np.arange(self.center[0]-math.floor(self.dims[0]/2),
                                 self.center[0]+math.floor(self.dims[0]/2)+1)
            y_indices = np.arange(self.center[1]-math.floor(self.dims[1]/2),
                                 self.center[1]+math.floor(self.dims[1]/2)+1)

        else:
            x_sums, y_sums = self.get_sums(img=False, array=array)
            x_indices = np.arange(len(x_sums))
            y_indices = np.arange(len(y_sums))

        x_centroid = 0
        y_centroid = 0
        total = 0
        for i in range(len(x_indices)):
            x_centroid = x_centroid + (x_indices[i]*x_sums[i])
            total = total + x_sums[i]
        for j in range(len(y_indices)):
            y_centroid = y_centroid + (y_indices[j]*y_sums[j])
        x_centroid = x_centroid / total
        y_centroid = y_centroid / total
        return x_centroid, y_centroid

    def get_uncertainty(self, img=True, array = []):
        if img:
            x_sums, y_sums = self.get_sums()
            x_centroid, y_centroid = self.get_centroids()
            x_indices = np.arange(self.center[0]-math.floor(self.dims[0]/2),
                                 self.center[0]+math.floor(self.dims[0]/2)+1)
            y_indices = np.arange(self.center[1]-math.floor(self.dims[1]/2),
                                 self.center[1]+math.floor(self.dims[1]/2)+1)
        else:
            x_sums, y_sums = self.get_sums(img=False, array = array)
            x_centroid, y_centroid = self.get_centroids(img=False, array=array)
            x_indices = np.arange(len(x_sums))
            y_indices = np.arange(len(y_sums))
        total = 0
        sig_x = 0
        sig_y = 0
        for i in range(len(x_indices)):
            sig_x = sig_x + (x_sums[i]*((x_indices[i]-x_centroid)**2))
            total = total + x_sums[i]
        for j in range(len(y_indices)):
            sig_y = sig_y + (y_sums[j]*((y_indices[j]-y_centroid)**2))
        sig_x = math.sqrt(sig_x/(total-1))
        sig_y = math.sqrt(sig_y/(total-1))
        sig_xbar = sig_x/math.sqrt(total)
        sig_ybar = sig_y/math.sqrt(total)
        return sig_xbar, sig_ybar



#centroid = Centroid('29p-data/20190901a/R/20190901a.R.0147',center=[1845,1530],dims=[50,50])
#centroid.display()
#print('x:' + str(centroid.get_centroids(img=True)[0]), 'y:' + str(centroid.get_centroids(img=True)[1]))
#print(centroid.get_uncertainty())
