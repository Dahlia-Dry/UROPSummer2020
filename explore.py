from starfield import StarField
from utilities import *
from astropy.io import fits
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from astropy.visualization import ZScaleInterval, LinearStretch,ImageNormalize
import numpy as np
import os
from itertools import compress
from mpl_toolkits.mplot3d import Axes3D
from centroid import Centroid
from scipy.optimize import curve_fit

#WEEK 1: QUALITATIVE EXPLORATION OF 29P DATA
#completed hand centroiding of 5 nights and checked results against mathematica
#and sextractor centroids
def mark_mm_centroids():
    centroids = pd.read_csv('29p-data/29P_Mma_Centers.tsv', sep='\t')
    targets = {}
    for i in range(len(centroids)):
        if centroids['file'].iloc[i][:9] not in targets:
            targets[centroids['file'].iloc[i][:9]] = []
        if centroids['file'].iloc[i][-4:] not in targets[centroids['file'].iloc[i][:9]]:
            targets[centroids['file'].iloc[i][:9]].append(centroids['file'].iloc[i][-4:])
    #print(targets)
    for date in targets:
        for number in targets[date]:
            data, fits_file, centroid = setuplmi(date,number,centroid=True,use_preset_path=True)
            find_obj(data, fits_file, 'centroid-offsets', centroid=centroid)

def mark_se_centroids(date):
    numbers = os.listdir('29p-data/' + date + '/' + 'R')
    bools = [i[0] != '.' and i[-4:] != 'fits' for i in numbers]
    numbers = list(compress(numbers, bools))
    numbers = [numbers[i][-4:] for i in range(len(numbers))]
    print(numbers)
    for number in numbers:
        findstar(date, number,destdir='29p-data/',designation='sextractor-offsets',
                use_preset_path=True)
#examples:
#mark_se_centroids('20191022a')
#_______________________________________________________________________________
#WEEK 2: First attempt at fitting + plotting Gaussian PSF
def subfield_select(image_data):
    def on_xlims_change(event_ax):
        global xlim
        xlim = event_ax.get_xlim()
        #print("updated xlims: ", event_ax.get_xlim())

    def on_ylims_change(event_ax):
        global ylim
        ylim = event_ax.get_ylim()
        #print("updated ylims: ", event_ax.get_ylim())

    fig = plt.figure()
    ax = fig.add_subplot(111)
    norm = ImageNormalize(image_data, interval=ZScaleInterval(),
                      stretch=LinearStretch())
    #ax.imshow(image_data, cmap = 'viridis',norm=norm)
    ax.imshow(image_data,vmin = image_data.mean(),
               vmax = 2*image_data.mean(),cmap='viridis')
    #ax.imshow(image_data, cmap = 'viridis')
    ax.callbacks.connect('xlim_changed', on_xlims_change)
    ax.callbacks.connect('ylim_changed', on_ylims_change)
    print('Close plot when selected window is satisfactory')
    plt.show()
    Xlim = [xlim[0],xlim[1]]
    Ylim = [ylim[0],ylim[1]]
    """if abs(int(Xlim[0])-int(Xlim[1])) % 2 != 0:
        Xlim[1] = Xlim[1]+1
    if abs(int(Ylim[0])-int(Ylim[1])) % 2 != 0:
        Ylim[1] = Ylim[1]+1"""
    print('X length:',str(abs(Xlim[0]-Xlim[1])))
    print('Y length:',str(abs(Ylim[0]-Ylim[1])))
    #Debug--
    image_data = image_data[int(Ylim[1]):int(Ylim[0]),int(Xlim[0]):int(Xlim[1])]
    norm = ImageNormalize(image_data, interval=ZScaleInterval(),
                      stretch=LinearStretch())
    print(image_data.shape)
    return image_data,Xlim,Ylim

def rebin(arr, new_shape):
    """Rebin 2D array arr to shape new_shape by averaging."""
    shape = (new_shape[0], arr.shape[0] // new_shape[0],
             new_shape[1], arr.shape[1] // new_shape[1])
    return arr.reshape(shape).mean(-1).mean(1)

def surface_plot(date,number,destdir,designation='contour-plots',
                plottype = 'surface',save=True):
    #select target object by manipulating zoom window
    #plot surface or contour of selected subframe
    rpath = os.path.join(destdir,date+'/R/'+date+'.R.'+number)
    rbpath = os.path.join(destdir,date+'/Rb/'+date+'.Rb.'+number)
    data = readcsv(rbpath)
    fits_file = rpath
    image_data = fits.getdata(fits_file)
    image_data,Xlim,Ylim = subfield_select(image_data)
    #plt.imshow(image_data, cmap = 'viridis',norm=norm)
    fig = plt.figure()
    fig, axes = plt.subplots(nrows=1, ncols=1, subplot_kw={'projection': '3d'})
    bnx, bny = (int(max(Xlim))-int(min(Xlim))), (int(max(Ylim))-int(min(Ylim)))
    x,y = np.linspace(1,bnx,bnx), np.linspace(1,bny,bny)
    x, y = np.meshgrid(x,y)
    z = rebin(image_data, (bny,bnx))
    print(x.shape,y.shape,z.shape)
    if plottype == 'surface':
        axes.plot_surface(x, y, z, rstride=4, cstride=4, linewidth=0.1,
                    antialiased=True, cmap=cm.plasma)
    elif plottype == 'contour':
        axes.contour(x, y, z, rstride=4, cstride=4, linewidth=0.1,
                    antialiased=True, cmap=cm.plasma)
    fig3 = plt.figure()
    ax = fig3.add_subplot(111)
    #ax.imshow(image_data, cmap = 'viridis',norm=norm)
    ax.imshow(image_data, cmap = 'viridis')
    if save:
        if not os.path.exists(str(os.path.join(destdir,date))):
            os.system('mkdir ' +str(os.path.join(destdir,date)))
            os.system('mkdir ' +str(os.path.join(destdir,date + '/'+designation)))
        elif not os.path.exists(str(os.path.join(destdir,date+'/'+designation))):
             os.system('mkdir ' +str(os.path.join(destdir,date + '/'+designation)))
        plt.savefig(str(os.path.join(destdir,date+'/'+designation+'/'+'contour'+
                        '_' + number+'.png')))
    plt.show()

def gaussian(coords,amplitude,xcentroid,ycentroid,sigma_x,sigma_y,offset):
    xcentroid = float(xcentroid)
    ycentroid = float(ycentroid)
    return offset+amplitude*np.exp(-(((coords[0]-xcentroid)**(2)/(2*sigma_x**(2)))+
                                    ((coords[1]-ycentroid)**(2)/(2*sigma_y**(2)))))

def fit_gaussian_v1(date,number,destdir,plottype=None,save=True,use_centroid=False):
    rpath = os.path.join(destdir,date+'/R/'+date+'.R.'+number)
    rbpath = os.path.join(destdir,date+'/Rb/'+date+'.Rb.'+number)
    data = readcsv(rbpath)
    fits_file = rpath
    image_data = fits.getdata(fits_file)
    image_data,Xlim,Ylim = subfield_select(image_data)
    norm = ImageNormalize(image_data, interval=ZScaleInterval(),
                      stretch=LinearStretch())
    maxindex = np.unravel_index(image_data.argmax(),image_data.shape)
    print(maxindex,image_data[maxindex[0]][maxindex[1]])
    if use_centroid:
        centroid = Centroid(rpath,center=[maxindex],dims=[20,20]).get_centroids()
        #plt.imshow(fits.getdata(fits_file),norm=norm)
        #plt.scatter(centroid[1]+Xlim[0],centroid[0]+Ylim[1],c='r')
        #plt.scatter(maxindex[1]+Xlim[0],maxindex[0]+Ylim[1],c='g')
        #plt.show()
        xcentroid = centroid[0]+Ylim[1]
        ycentroid = centroid[1]+Xlim[0]
    else:
        xcentroid = maxindex[0]+Ylim[1]
        ycentroid = maxindex[1]+Xlim[0]
    print('x:',xcentroid,'y:',ycentroid)
    bnx, bny = (int(max(Xlim))-int(min(Xlim))), (int(max(Ylim))-int(min(Ylim)))
    x,y = np.linspace(1,bnx,bnx), np.linspace(1,bny,bny)
    x, y = np.meshgrid(x,y)
    z = image_data
    #print(x.shape,y.shape, image_data.shape)
    #z = rebin(image_data, (bny,bnx))
    #axes.plot_wireframe(x,y,z)
    #plt.show()
    xflat = x.flatten()
    yflat = y.flatten()
    zflat = z.flatten()
    popt, pcov = curve_fit(gaussian, (xflat,yflat),zflat)
    print(popt)
    data = [x,y,z]
    fit = [popt, pcov]
    if plottype is not None:
        plot_v1fit(plottype,data,fit,date,number,destdir,save=save)
    return data, fit

def plot_v1fit(plottype,data,fit,date,number,destdir,
                designation='gaussian-fits-v1',save=True):
    x = data[0]
    y = data[1]
    z = data[2]
    popt = fit[0]
    pcov = fit[1]
    xflat = x.flatten()
    yflat = y.flatten()
    zflat = z.flatten()
    if plottype is '3dmesh':
        fig, axes = plt.subplots(nrows=1, ncols=1, subplot_kw={'projection': '3d'})
        axes.set_title('2D Gaussian PSF: ' +
            r'$G(x,y) = \beta + \alpha*e^{-(\frac{(x-x_0)^2}{2*\sigma_x}+\frac{(y-y_0)^2}{2*\sigma_y})}$')
        params = '\n'.join((
            r'$\alpha=%.2f$' % (popt[0], ),
            r'$\beta=%.2f$' % (popt[-1], ),
            r'$x_0=%.2f$' % (popt[1], ),
            r'$y_0=%.2f$' % (popt[2], ),
            r'$\sigma_x=%.2f$' % (popt[3], ),
            r'$\sigma_y=%.2f$' % (popt[4], ),))
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        axes.text2D(0.05, 0.05, params, transform=axes.transAxes, fontsize=8,
                    verticalalignment='bottom', bbox=props)
        axes.plot_wireframe(x,y,np.reshape(gaussian((xflat,yflat),*popt),(x.shape)))
        axes.contour(x, y, z, rstride=4, cstride=4, linewidth=0.1,
                    antialiased=True,cmap=cm.plasma)
    elif plottype is '3dcontour':
        fig, axes = plt.subplots(nrows=1, ncols=1, subplot_kw={'projection': '3d'})
        axes.set_title('2D Gaussian PSF: ' +
            r'$G(x,y) = \beta + \alpha*e^{-(\frac{(x-x_0)^2}{2*\sigma_x}+\frac{(y-y_0)^2}{2*\sigma_y})}$')
        params = '\n'.join((
            r'$\alpha=%.2f$' % (popt[0], ),
            r'$\beta=%.2f$' % (popt[-1], ),
            r'$x_0=%.2f$' % (popt[1], ),
            r'$y_0=%.2f$' % (popt[2], ),
            r'$\sigma_x=%.2f$' % (popt[3], ),
            r'$\sigma_y=%.2f$' % (popt[4], ),))
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        axes.text2D(0.05, 0.05, params, transform=axes.transAxes, fontsize=8,
                    verticalalignment='bottom', bbox=props)
        axes.contour(x,y,np.reshape(gaussian((xflat,yflat),*popt),(x.shape)),
                    rstride=4, cstride=4, linewidth=0.1,antialiased=True,
                    cmap=cm.viridis)
        axes.contour(x, y, z, rstride=4, cstride=4, linewidth=0.1,
                    antialiased=True,cmap=cm.plasma)
    elif plottype is '2d':
        fig,ax = plt.subplots()
        ax.set_title('2D Gaussian PSF: ' +
            r'$G(x,y) = \beta + \alpha*e^{-(\frac{(x-x_0)^2}{2*\sigma_x}+\frac{(y-y_0)^2}{2*\sigma_y})}$')
        params = '\n'.join((
            r'$\alpha=%.2f$' % (popt[0], ),
            r'$\beta=%.2f$' % (popt[-1], ),
            r'$x_0=%.2f$' % (popt[1], ),
            r'$y_0=%.2f$' % (popt[2], ),
            r'$\sigma_x=%.2f$' % (popt[3], ),
            r'$\sigma_y=%.2f$' % (popt[4], ),))
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        ax.text(0.05, 0.05, params, transform=ax.transAxes, fontsize=8,
                    verticalalignment='bottom', bbox=props)
        rawdata = ax.contour(x,y,z, cmap=cm.Blues)
        fitteddata = ax.contour(x,y,np.reshape(gaussian((xflat,yflat),*popt),(x.shape)),
                                cmap=cm.Greens)
        plt.legend([rawdata,fitteddata],['Raw Data','Fitted Curve'])
    elif plottype is 'resid':
        fig,ax = plt.subplots()
        ax.set_title('Residuals for 2D Gaussian PSF: ' +
            r'$G(x,y) = \beta + \alpha*e^{-(\frac{(x-x_0)^2}{2*\sigma_x}+\frac{(y-y_0)^2}{2*\sigma_y})}$')
        params = '\n'.join((
            r'$\alpha=%.2f$' % (popt[0], ),
            r'$\beta=%.2f$' % (popt[-1], ),
            r'$x_0=%.2f$' % (popt[1], ),
            r'$y_0=%.2f$' % (popt[2], ),
            r'$\sigma_x=%.2f$' % (popt[3], ),
            r'$\sigma_y=%.2f$' % (popt[4], ),))
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        ax.text(0.05, 0.05, params, transform=ax.transAxes, fontsize=8,
                    verticalalignment='bottom', bbox=props)
        resids = z - np.reshape(gaussian((xflat,yflat),*popt),(x.shape))
        print(np.arange(len(resids)).shape,resids.shape)
        plt.imshow(resids)
    else:
        raise Exception("Please Specify plot type: 3dmesh, 3dcontour,\
                        2d, resid")

    if save:
        if not os.path.exists(str(os.path.join(destdir,date))):
            os.system('mkdir ' +str(os.path.join(destdir,date)))
            os.system('mkdir ' +str(os.path.join(destdir,date + '/'+designation)))
        elif not os.path.exists(str(os.path.join(destdir,date+'/'+designation))):
             os.system('mkdir ' +str(os.path.join(destdir,date + '/'+designation)))
        plt.savefig(str(os.path.join(destdir,date+'/'+designation+'/'+'gaussianfit'+
                        '-' + plottype + '_' +number+'.png')))
    plt.show()

#surface_plot('20190901a','0147',destdir='29p-data/',save=False)
#data, fit = fit_gaussian_v1('20190901a','0147',destdir='29p-data/',plottype='resid')
#___________________________________________________________________________________
#WEEK 3:
def fit_gaussian_v2(date,number,destdir,plottype=None,save=True,use_centroid=False):
    rpath = os.path.join(destdir,date+'/R/'+date+'.R.'+number)
    rbpath = os.path.join(destdir,date+'/Rb/'+date+'.Rb.'+number)
    data = readcsv(rbpath)
    fits_file = rpath
    image_data = fits.getdata(fits_file)
    image_data,Xlim,Ylim = subfield_select(image_data)
    norm = ImageNormalize(image_data, interval=ZScaleInterval(),
                      stretch=LinearStretch())
    maxindex = np.unravel_index(image_data.argmax(),image_data.shape)
    
