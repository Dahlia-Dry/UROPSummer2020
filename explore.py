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
import matplotlib.patches as patches

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

def fit_gaussian_v1(date,number,destdir,plottype=None, subfield=None,lims=None,
                    save=True,use_centroid=False):
    rpath = os.path.join(destdir,date+'/R/'+date+'.R.'+number)
    rbpath = os.path.join(destdir,date+'/Rb/'+date+'.Rb.'+number)
    data = readcsv(rbpath)
    fits_file = rpath
    image_data = fits.getdata(fits_file)
    if subfield is None:
        image_data,Xlim,Ylim = subfield_select(image_data)
    else:
        image_data = subfield
        Xlim, Ylim = lims[0],lims[1]
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
        if plottype is not 'all':
            plot_v1fit(plottype,data,fit,date,number,destdir,save=save)
        else:
            fig = plt.figure()
            fig.suptitle('Gaussian PSF fit for '+date+'-'+number+':'+
            r'$G(x,y) = \beta + \alpha*e^{-(\frac{(x-x_0)^2}{2*\sigma_x}+\frac{(y-y_0)^2}{2*\sigma_y})}$')
            ax = fig.add_subplot(221)
            ax.set_title('Object')
            ax.axes.xaxis.set_visible(False)
            ax.axes.yaxis.set_ticks([])
            ax.imshow(image_data, vmin = image_data.mean(),
                       vmax = 2*image_data.mean(),cmap='viridis')
            ax = fig.add_subplot(222,projection='3d')
            ax.set_title('Model Fit')
            plot_v1fit('3dmesh',data,fit,date,number,destdir,designation='gaussian-fits-v1',
                        figax = [fig,ax],showparams=False,save=save)
            ax = fig.add_subplot(223)
            ax.axes.xaxis.set_visible(False)
            ax.axes.yaxis.set_ticks([])
            ax.set_title('Observed vs Model Contour')
            plot_v1fit('2d',data,fit,date,number,destdir,designation='gaussian-fits-v1',
                        figax = [fig,ax],smalltext=True,showparams=True,save=save)
            ax = fig.add_subplot(224)
            ax.axes.xaxis.set_ticks([])
            ax.axes.yaxis.set_ticks([])
            ax.set_title('Residuals')
            plot_v1fit('resid',data,fit,date,number,destdir,designation='gaussian-fits-v1',
                        figax = [fig,ax],showparams=False,save=save)
            plt.show()
    return data, fit

def plot_v1fit(plottype,data,fit,date,number,destdir,
                designation='gaussian-fits-v1',showparams=True,
                figax = None,smalltext=False, save=True):
    x = data[0]
    y = data[1]
    z = data[2]
    popt = fit[0]
    pcov = fit[1]
    xflat = x.flatten()
    yflat = y.flatten()
    zflat = z.flatten()
    if plottype is '3dmesh':
        if figax is None:
            fig, ax = plt.subplots(nrows=1, ncols=1, subplot_kw={'projection': '3d'})
            ax.set_title('2D Gaussian PSF: ' +
                r'$G(x,y) = \beta + \alpha*e^{-(\frac{(x-x_0)^2}{2*\sigma_x}+\frac{(y-y_0)^2}{2*\sigma_y})}$')
        else:
            fig, ax = figax[0], figax[1]
        if showparams:
            params = '\n'.join((
                r'$\alpha=%.2f$' % (popt[0], ),
                r'$\beta=%.2f$' % (popt[-1], ),
                r'$x_0=%.2f$' % (popt[1], ),
                r'$y_0=%.2f$' % (popt[2], ),
                r'$\sigma_x=%.2f$' % (popt[3], ),
                r'$\sigma_y=%.2f$' % (popt[4], ),))
            props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
            if smalltext:
                fontsize=4
            else:
                fontsize=8
            ax.text2D(0.05, 0.05, params, transform=ax.transAxes, fontsize=fontsize,
                        verticalalignment='bottom', bbox=props)
        ax.plot_wireframe(x,y,np.reshape(gaussian((xflat,yflat),*popt),(x.shape)))
        ax.contour(x, y, z, rstride=4, cstride=4, linewidth=0.1,
                    antialiased=True,cmap=cm.plasma)
    elif plottype is '3dcontour':
        if figax is None:
            fig, ax = plt.subplots(nrows=1, ncols=1, subplot_kw={'projection': '3d'})
            ax.set_title('2D Gaussian PSF: ' +
                r'$G(x,y) = \beta + \alpha*e^{-(\frac{(x-x_0)^2}{2*\sigma_x}+\frac{(y-y_0)^2}{2*\sigma_y})}$')
        else:
            fig, ax = figax[0], figax[1]
        if showparams:
            params = '\n'.join((
                r'$\alpha=%.2f$' % (popt[0], ),
                r'$\beta=%.2f$' % (popt[-1], ),
                r'$x_0=%.2f$' % (popt[1], ),
                r'$y_0=%.2f$' % (popt[2], ),
                r'$\sigma_x=%.2f$' % (popt[3], ),
                r'$\sigma_y=%.2f$' % (popt[4], ),))
            props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
            if smalltext:
                fontsize=4
            else:
                fontsize=8
            ax.text2D(0.05, 0.05, params, transform=ax.transAxes, fontsize=fontsize,
                        verticalalignment='bottom', bbox=props)
        ax.contour(x,y,np.reshape(gaussian((xflat,yflat),*popt),(x.shape)),
                    rstride=4, cstride=4, linewidth=0.1,antialiased=True,
                    cmap=cm.viridis)
        ax.contour(x, y, z, rstride=4, cstride=4, linewidth=0.1,
                    antialiased=True,cmap=cm.plasma)
    elif plottype is '2d':
        if figax is None:
            fig, ax = plt.subplots()
            ax.set_title('2D Gaussian PSF: ' +
                r'$G(x,y) = \beta + \alpha*e^{-(\frac{(x-x_0)^2}{2*\sigma_x}+\frac{(y-y_0)^2}{2*\sigma_y})}$')
        else:
            fig, ax = figax[0], figax[1]
        if showparams:
            params = '\n'.join((
                r'$\alpha=%.2f$' % (popt[0], ),
                r'$\beta=%.2f$' % (popt[-1], ),
                r'$x_0=%.2f$' % (popt[1], ),
                r'$y_0=%.2f$' % (popt[2], ),
                r'$\sigma_x=%.2f$' % (popt[3], ),
                r'$\sigma_y=%.2f$' % (popt[4], ),))
            props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
            if smalltext:
                fontsize=4
            else:
                fontsize=8
            ax.text(0.05, 0.05, params, transform=ax.transAxes, fontsize=fontsize,
                        verticalalignment='bottom', bbox=props)
        rawdata = ax.contour(x,y,z, cmap=cm.Blues)
        fitteddata = ax.contour(x,y,np.reshape(gaussian((xflat,yflat),*popt),(x.shape)),
                                cmap=cm.Greens)
        #plt.legend([rawdata,fitteddata],['Raw Data','Fitted Curve'])
    elif plottype is 'resid':
        if figax is None:
            fig, ax = plt.subplots()
            ax.set_title('Residuals for 2D Gaussian PSF: ' +
                r'$G(x,y) = \beta + \alpha*e^{-(\frac{(x-x_0)^2}{2*\sigma_x}+\frac{(y-y_0)^2}{2*\sigma_y})}$')
        else:
            fig, ax = figax[0], figax[1]
        if showparams:
            params = '\n'.join((
                r'$\alpha=%.2f$' % (popt[0], ),
                r'$\beta=%.2f$' % (popt[-1], ),
                r'$x_0=%.2f$' % (popt[1], ),
                r'$y_0=%.2f$' % (popt[2], ),
                r'$\sigma_x=%.2f$' % (popt[3], ),
                r'$\sigma_y=%.2f$' % (popt[4], ),))
            props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
            if smalltext:
                fontsize=4
            else:
                fontsize=8
            ax.text(0.05, 0.05, params, transform=ax.transAxes, fontsize=fontsize,
                        verticalalignment='bottom', bbox=props)
        resids = z - np.reshape(gaussian((xflat,yflat),*popt),(x.shape))
        print(np.arange(len(resids)).shape,resids.shape)
        ax.imshow(resids,cmap=cm.plasma)
        color = ax.imshow(resids,cmap=cm.plasma)
        fig.colorbar(color)
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

    if figax is not None:
        pass
    else:
        plt.show()


#surface_plot('20190901a','0147',destdir='29p-data/',save=False)
#data, fit = fit_gaussian_v1('20190901a','0147',destdir='29p-data/',plottype='all')
#___________________________________________________________________________________
#WEEK 3:
def fit_gaussian_v2(date,number,destdir,starap = 15, objap = 20,
                    plottype=None,save=True,use_centroid=False):
    rpath = os.path.join(destdir,date+'/R/'+date+'.R.'+number)
    rbpath = os.path.join(destdir,date+'/Rb/'+date+'.Rb.'+number)
    data = readcsv(rbpath)
    fits_file = rpath
    image_data = fits.getdata(fits_file)
    #step 1: choose subfield
    print('Select a window that includes the object and at least 3 calibration stars:')
    image_data,Xlim,Ylim = subfield_select(image_data)
    #print('Xlim:',Xlim)
    #print('Ylim:',Ylim)
    #step 2: pick stars to calculate sigma with
    print('Pick out 3 calibration stars;\
            window will close automatically once selections are made')
    rbindices = findstar(rpath=rpath,rbpath = rbpath, destdir = '29p-data/',
                        nstars = 3, subfield = image_data,lims=[Xlim,Ylim], save=False)
    #step 3: isolate obj in its own subfield
    print('Use zoom to select a subfield containing only the object of interest;\
            close window when done')
    obj_data,objX, objY = subfield_select(image_data)
    fig = plt.figure()
    #Subplot 1 of 4: show apertures of fit windows
    ax = fig.add_subplot(2,2,1)
    ax.imshow(image_data, vmin = image_data.mean(),
               vmax = 2*image_data.mean(),cmap='viridis')
    slx,sly,subfields = [], [], []
    for index in rbindices:
        #print('rb:',data.iloc[index]['XIM'],data.iloc[index]['YIM'])
        coords = (data.iloc[index]['XIM']-Xlim[0]-starap/2,
                    data.iloc[index]['YIM']-Ylim[1]-starap/2)
        slx.append([coords[0]-starap/2,coords[0]+starap/2])
        sly.append([coords[1]+starap/2,coords[1]-starap/2])
        subfields.append(image_data[int(coords[1]):int(coords[1]+starap),
                                    int(coords[0]):int(coords[0]+starap)])
        #print('coords:',coords)
        rect = patches.Rectangle(coords,starap, starap, edgecolor='r', fill=False)
        ax.add_patch(rect)
    #print(np.amax(obj_data),np.amax(image_data))
    objcoords = list(zip(*np.where(image_data == np.amax(obj_data))))
    objcoords = (objcoords[0][1]-objap/2,
                objcoords[0][0]-objap/2)
    obx = [objcoords[0]-objap/2,objcoords[0]+objap/2]
    oby = [objcoords[1]+objap/2,objcoords[1]-objap/2]
    objfield = image_data[int(objcoords[1]):int(objcoords[1]+objap),
                                int(objcoords[0]):int(objcoords[0]+objap)]
    print(objfield.size)
    rect = patches.Rectangle(objcoords, objap, objap, edgecolor = 'black',fill=False)
    ax.add_patch(rect)
    ax.set_title('Fit Windows')
    ax = fig.add_subplot(2,2,2)
    ax.imshow(objfield)
    #get calibration sigmas
    data, fit = [],[]
    fig2 = plt.figure()
    fig2.suptitle('Calibration Star Fits:'+r'$G(x,y) = \beta + \alpha*e^{-(\frac{(x-x_0)^2}{2*\sigma_x}+\frac{(y-y_0)^2}{2*\sigma_y})}$')
    for i in range(len(subfields)):
        d, f = fit_gaussian_v1(date,number,destdir=destdir,
                                    subfield=subfields[i],lims=[slx[i],sly[i]])
        data.append(d)
        fit.append(f)
        #plot subfields
        sax = fig2.add_subplot(3,3,i+1)
        sax.axes.xaxis.set_visible(False)
        sax.axes.yaxis.set_ticks([])
        sax.imshow(subfields[i])
        sax.set_title('Star ' + str(i+1))
        #plot residuals
        resid =fig2.add_subplot(3,3,i+4)
        resid.axes.xaxis.set_visible(False)
        resid.axes.yaxis.set_ticks([])
        plot_v1fit('resid',d,f,date,number,destdir,
                        designation='gaussian-fits-v1',
                        figax = [fig2,resid],
                        showparams=False,save=False)
        #plot contours
        contour = fig2.add_subplot(3,3,i+7)
        contour.axes.xaxis.set_visible(False)
        contour.axes.yaxis.set_ticks([])
        plot_v1fit('2d',d,f,date,number,destdir,
                        designation='gaussian-fits-v1',
                        figax = [fig2,contour],
                        showparams=True,smalltext=True,save=False)
        if i == 0:
            sax.set_ylabel('Subfield')
            resid.set_ylabel('Residuals')
            contour.set_ylabel('Contours')
    plt.show()
    sigmas = [fit[i][0][4] for i in range(len(subfields))]
    sigma_guess= sum(sigmas)/len(sigmas)
    beta_guess = np.median(objfield)

fit_gaussian_v2('20190901a','0147',destdir='29p-data/')
