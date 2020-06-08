from starfield import StarField
from utilities import *
from astropy.io import fits
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from astropy.visualization import ZScaleInterval, LinearStretch,ImageNormalize
import numpy as np
import os

def find_obj(data, image_data, centroid=None, destdir = '29p-data/',
            save=True, size=100):
    def onclick(event):
        global ix, iy
        ix, iy = event.xdata, event.ydata
        print('x = %d, y = %d'%(ix, iy))
        fig.canvas.mpl_disconnect(cid)

    def find_nearest(array, value):
        idx = (np.abs(array-value)).argmin()
        return idx
    norm = ImageNormalize(image_data, interval=ZScaleInterval(),
                      stretch=LinearStretch())
    #Allow user to pick the centroid of their object
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.imshow(image_data, cmap = 'viridis',norm=norm)
    #fig.colorbar()
    cid = fig.canvas.mpl_connect('button_press_event', onclick)
    plt.show()

    if centroid is None:
        #Find closest entry in .Rb file to user click by
        #looking for where distance is closest to 0
        distances = np.sqrt((data['XIM'].values-ix)**2 +(data['YIM'].values-iy)**2)
        closest = find_nearest(distances, 0)
        centroid = [data['XIM'].iloc[closest], data['YIM'].iloc[closest]]
        distance=distances[closest]
    print('centroid:' + str(centroid))
    distance= np.sqrt((ix-centroid[0])**2+(iy-centroid[1])**2)
    plt.close(fig)
    fig2 = plt.figure()
    ax2 = fig2.add_subplot(111)
    bounds = [int(iy-size),int(iy+size),int(ix-size),int(ix+size)]
    ax2.imshow(image_data,cmap = 'viridis',norm=norm)
    ax2.plot([ix,centroid[0]],[iy,centroid[1]],c='r')
    if save:
        if not os.path.exists(destdir+date):
            os.system('mkdir ' +destdir+date)
            os.system('mkdir ' +destdir+date + '/centroid-offsets')
        elif not os.path.exists(destdir+date+'/centroid-offsets'):
             os.system('mkdir ' +destdir+date + '/centroid-offsets')
        plt.savefig(destdir+date+'/centroid-offsets/'+'offset='+
                    str(int(distance))+ '_' + number+'.png')
    print(data[closest])
    plt.show()

def mark_centroids():
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
            data, fits_file, centroid = setup(date,number,centroid=True,use_preset_path=True)
            find_obj(data, fits_file, centroid=centroid)
