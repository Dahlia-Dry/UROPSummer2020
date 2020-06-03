from astropy.io import fits
import pandas as pd
from utilities import readcsv
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os.path
import numpy as np
import scipy.stats
import pickle
import ipyvolume as p3
import matplotlib.cm as cm
import matplotlib as mat
import ipywidgets as widgets
from IPython.display import display
from ipywidgets import FloatSlider, ColorPicker, VBox, jslink
from scipy import stats


class StarField(object):
    """
    Object representing an individual starfield in the mitastrometry database.
    Defined by the following arguments:
    -telescope (str) = telescope used, eg. 'SARAS'
    -date (str) = date taken, eg. '20190823'
    -number (str) = frame number, eg. '0357'
    -data (Pandas DataFrame) = .Rb data for each star in image, can be read using csv_reader helper function
    -fits_file (str) = name of .R file, does not need .fits extension
    -create (bool) = True if you want to make a corresponding directory for the starfield
    """
    def __init__(self, telescope, date, number, data, fits_file, create=True):
        self.telescope = telescope
        self.date = date
        self.number = number
        self.data = data
        self.create= create
        self.fits_file = fits_file
        self.directory = 'mitastrometry/' + str(self.telescope) + '/' + str(self.date) + '/' + str(self.date) + '_' + str(self.number) + "_analysis"
        """else:
            raise Exception('Directory for ' + self.date + '_' + self.number +
                            "already exists. If you wish to perform analysis again, delete existing directory.")"""
        self.outfile = os.path.join(self.directory, self.number + 'analysis.txt')
        #self.instrument = fits.open(fits_file)[0].header['INSTRUME']
        self.object = fits.open(fits_file)[0].header['OBJECT']

    def describe(self):
        return {'telescope':self.telescope,'date': self.date,'number': self.number,
                'object':self.object}

    def write_outfile(self, mpv=False):
        if self.create:
            if not os.path.exists(self.directory):
                os.makedirs(self.directory)
        if mpv:
            if not os.path.exists(os.path.join(self.directory, self.number + '_mpvs.csv')):
                mpv_list = pd.Series([self.mpv_excluding_annulus(star) for star in range(len(self.data['Eta']))],
                                     index=None)
                mpv_list.to_csv(os.path.join(self.directory, self.number + '_mpvs.csv'))
            else:
                mpv_list = pd.read_csv(os.path.join(self.directory, self.number + '_mpvs.csv'))
        n = 1
        f = open(self.outfile, 'w')
        f.write('____________________________________________________________________\n')
        print('____________________________________________________________________\n')
        f.write('STATISTICAL ANALYSIS OF ' + self.date + "_" + self.number + "\n")
        print('STATISTICAL ANALYSIS OF ' + self.date + "_" + self.number + "\n")
        f.write('telescope: ' + self.telescope + '\n')
        print('telescope: ' + self.telescope + '\n')
        #f.write('instrument: ' + str(self.instrument) + '\n')
        #print('instrument: ' + str(self.instrument) + '\n')
        f.write('____________________________________________________________________\n')
        print('____________________________________________________________________\n')
        f.write('RESIDUAL PLOTS\n')
        print('RESIDUAL PLOTS\n')

        f.write('Xi Residual Magnitudes:\n')
        print('Xi Residual Magnitudes:\n')
        ximu, xisig, r, r2 = self.plot_resids('xi', 'fig'+str(n))
        f.write('fig' + str(n) + '\t' + 'mean = ' + str(ximu) + '\t std = ' + str(xisig) + '\n')
        print('fig' + str(n) +'\t' + 'mean = ' + str(ximu) + '\t std = ' + str(xisig) + '\n')
        f.write('r = ' + str(r) + '\t r^2 = ' + str(r2) + '\n')
        print('r = ' + str(r) + '\t r^2 = ' + str(r2) + '\n')
        n = n + 1

        f.write('Eta Residual Magnitudes:\n')
        print('Eta Residual Magnitudes:\n')
        etamu, etasig, r, r2 = self.plot_resids('eta', 'fig' + str(n))
        f.write('fig' + str(n) +  '\t' + 'mean = ' + str(etamu) + '\t std = ' + str(etasig) + '\n')
        print('fig'  + str(n) + '\t' + 'mean = ' + str(etamu) + '\t std = ' + str(etasig) + '\n')
        f.write('r = ' + str(r) + '\t r^2 = ' + str(r2) + '\n')
        print('r = ' + str(r) + '\t r^2 = ' + str(r2) + '\n')
        n = n + 1

        f.write('____________________________________________________________________\n')
        print('____________________________________________________________________\n')
        f.write('UNIVARIATE ANALYSIS\n')
        print('UNIVARIATE ANALYSIS\n')

        f.write('Residual Magnitude:')
        print('Residual Magnitude:')
        mu, sig = self.univariate('resids', 'fig' + str(n))
        f.write('fig' + str(n) + '\t' + 'mean = ' + str(mu) + '\t std = ' + str(sig) + '\n')
        print('fig' + str(n) +  '\t' + 'mean = ' + str(mu) + '\t std = ' + str(sig) + '\n')
        n = n + 1

        if mpv:
            f.write('Subplot Mean Pixel Value:')
            print('Subplot Mean Pixel Value:')
            mu, sig = self.univariate('mpv', 'fig' + str(n))
            f.write('fig' + str(n) + '\t' + 'mean = ' + str(mu) + '\t std = ' + str(sig) + '\n')
            print('fig' + str(n) +'\t' + 'mean = ' + str(mu) + '\t std = ' + str(sig) + '\n')
            n = n + 1

        f.write('Stellar Magnitude:')
        print('Stellar Magnitude:')
        mu, sig = self.univariate('MagAuto', 'fig' + str(n))
        f.write('fig'  + str(n) + '\t' + 'mean = ' + str(mu) + '\t std = ' + str(sig) + '\n')
        print('fig' + str(n) + '\t' + 'mean = ' + str(mu) + '\t std = ' + str(sig) + '\n')
        n = n + 1

        f.write('FWHM:')
        print('FWHM:')
        mu, sig = self.univariate('FWHM', 'fig' + str(n))
        f.write('fig' + str(n) +  '\t' + 'mean = ' + str(mu) + '\t std = ' + str(sig) + '\n')
        print('fig' + str(n) + '\t' + 'mean = ' + str(mu) + '\t std = ' + str(sig) + '\n')
        n = n + 1

        f.write('____________________________________________________________________\n')
        print('____________________________________________________________________\n')
        f.write('BIVARIATE ANALYSIS\n')
        print('BIVARIATE ANALYSIS\n')
        f.write('--Colormap overlay is a Gaussian Kernel Density Estimation--\n')
        print('--Colormap overlay is a Gaussian Kernel Density Estimation--\n')
        f.write('--KDEs are zoomed to middle 80% of data for visual clarity--\n')
        print('--KDEs are zoomed to middle 80% of data for visual clarity--\n')
        f.write('--Raw scatterplots show totality of data--\n')
        print('--Raw scatterplots show totality of data--\n')

        f.write('Residual Magnitude vs Stellar Magnitude:\n')
        print('Residual Magnitude vs Stellar Magnitude:\n')
        f.write('figs ' + str(n) + ',' + str(n+1) + ':\n')
        print('figs ' + str(n) + ',' + str(n+1) + ':\n')
        r, r2 = self.bivariate('MagAuto','resids', 'fig' + str(n), 'fig' + str(n+1))
        f.write('r = ' + str(r) + '\t r^2 = ' + str(r2) + '\n')
        print('r = ' + str(r) + '\t r^2 = ' + str(r2) + '\n')
        n = n+2

        f.write('Residual Magnitude vs Subplot Mean Pixel Value:\n')
        print('Residual Magnitude vs Subplot Mean Pixel Value:\n')
        f.write('figs ' + str(n) + ',' + str(n+1) + ':\n')
        print('figs ' + str(n) + ',' + str(n+1) + ':\n')
        r, r2 = self.bivariate('mpv','resids', 'fig' + str(n), 'fig' + str(n+1))
        f.write('r = ' + str(r) + '\t r^2 = ' + str(r2) + '\n')
        print('r = ' + str(r) + '\t r^2 = ' + str(r2) + '\n')
        n = n + 2

        f.write('Residual Magnitude vs FWHM:\n')
        print('Residual Magnitude vs FWHM:\n')
        f.write('figs ' + str(n) + ',' + str(n+1) + ':\n')
        print('figs ' + str(n) + ',' + str(n+1) + ':\n')
        r, r2 = self.bivariate('FWHM', 'resids' ,'fig' + str(n), 'fig' + str(n+1))
        f.write('r = ' + str(r) + '\t r^2 = ' + str(r2) + '\n')
        print('r = ' + str(r) + '\t r^2 = ' + str(r2) + '\n')
        n = n + 2

        f.write('Stellar Magnitude vs FWHM:\n')
        print('Stellar Magnitude vs FWHM:\n')
        f.write('figs ' + str(n) + ',' + str(n + 1) + ':\n')
        print('figs ' + str(n) + ',' + str(n + 1) + ':\n')
        r, r2 = self.bivariate('MagAuto', 'FWHM', 'fig' + str(n), 'fig' + str(n + 1),diff=True)
        f.write('r = ' + str(r) + '\t r^2 = ' + str(r2) + '\n')
        print('r = ' + str(r) + '\t r^2 = ' + str(r2) + '\n')
        n = n + 2

        f.write('Elongation vs Ellipticity:\n')
        print('Elongation vs Ellipticity:\n')
        f.write('figs ' + str(n) + ',' + str(n + 1) + ':\n')
        print('figs ' + str(n) + ',' + str(n + 1) + ':\n')
        r, r2 = self.bivariate('Elong', 'Ellip', 'fig' + str(n), 'fig' + str(n + 1),diff=True)
        f.write('r = ' + str(r) + '\t r^2 = ' + str(r2) + '\n')
        print('r = ' + str(r) + '\t r^2 = ' + str(r2) + '\n')
        n = n + 2

        f.write('____________________________________________________________________\n')
        print('____________________________________________________________________\n')
        f.write('MULTIVARIATE ANALYSIS\n')
        print('MULTIVARIATE ANALYSIS\n')

        f.write('3D Residual Quiver Plot: fig' + str(n) +  '\n')
        print('3D Residual Quiver Plot: fig' + str(n) + '\n')
        self.plot_3D('quiver', 'fig' + str(n))
        n = n + 1

        f.write('3D Residual Surface Triangulation: fig' + str(n) + '\n')
        print('3D Residual Surface Triangulation: fig' + str(n) + '\n')
        self.plot_3D('trisurf', 'fig' + str(n))
        n = n + 1
        f.close()

    def getdata(self, param):
        if param == "resids":
            oceta = self.data['O-Ceta']
            oxi = self.data['O-Cxi']
            numbers = np.array([(((oceta[i] ** 2) + (oxi[i] ** 2)) ** .5) for i in range(len(oceta))])
        elif param == "mpv":
            numbers = np.array(pd.read_csv(os.path.join(self.directory, self.number + '_mpvs.csv'), header=None)[1])
        else:
            numbers = self.data[param]
        return numbers

    def plot_resids(self, param, fname):
        if param == 'xi':
            oxi = self.data['O-Cxi']
            xi = self.data['Xi']
            ximu = np.mean(oxi)
            xisig = np.std(oxi)
            r = np.corrcoef(oxi, xi)[0][1]
            plt.xlabel('Xi')
            plt.ylabel('O-Cxi')
            plt.title('Xi Residuals')
            plt.scatter(xi, oxi, s=1)
            plt.savefig(os.path.join(self.directory, fname + '.png'))
            plt.show()
            return ximu, xisig, r, r**2
        elif param == 'eta':
            oeta = self.data['O-Ceta']
            eta = self.data['Eta']
            etamu = np.mean(oeta)
            etasig = np.std(oeta)
            r= np.corrcoef(oeta, eta)[0][1]
            plt.xlabel('Eta')
            plt.ylabel('O-Ceta')
            plt.title('Eta Residuals')
            plt.scatter(eta, oeta, s=1, color = 'orange')
            plt.savefig(os.path.join(self.directory, fname + '.png'))
            plt.show()
            return etamu, etasig, r, r**2

    def univariate(self, param, fname):
        label = param
        if param == "resids":
            oceta = self.data['O-Ceta']
            oxi = self.data['O-Cxi']
            numbers = np.array([(((oceta[i] ** 2) + (oxi[i] ** 2)) ** .5) for i in range(len(oceta))])
            label = 'Residual Magnitude'
        elif param == "mpv":
            numbers = np.array(pd.read_csv(os.path.join(self.directory, self.number + '_mpvs.csv'), header=None)[1])
            label = 'Subplot Mean Pixel Value'
        else:
            numbers = self.data[param]

        bins = int(1+3.22*np.log(len(numbers))) #Sturge's Rule
        mu = np.mean(numbers)
        sig = np.std(numbers)
        plt.hist(numbers, bins)
        plt.title(label + " Distribution")
        plt.savefig(os.path.join(self.directory, fname + '.png'))
        plt.show()
        return mu, sig

    def bivariate(self, x,y, fname, fname2, diff=False):
        xlabel = x
        ylabel = y
        if y == "resids":
            oceta = self.data['O-Ceta']
            oxi = self.data['O-Cxi']
            param2_i = np.array([(((oceta[i] ** 2) + (oxi[i] ** 2)) ** .5) for i in range(len(oceta))])
            ylabel = 'Residual Magnitude'
        else:
            param2_i = np.array(self.data[y])

        if x == "mpv":
            param1_i = np.array(pd.read_csv(os.path.join(self.directory, self.number + '_mpvs.csv'), header=None)[1])
            #param1_i = [self.get_mean_subplot_value(i) for i in range(len(self.data['Eta']))]
            xlabel = 'Subplot Mean Pixel Value'
        else:
            param1_i = np.array(self.data[x])

        param1 = []
        param2 = []
        x1 = np.percentile(param1_i, 10)
        x3 = np.percentile(param1_i, 90)
        y1 = np.percentile(param2_i, 10)
        y3 = np.percentile(param2_i, 90)
        for i in range(len(param1_i)):
            if param1_i[i] > x1 and param1_i[i] < x3 and param2_i[i] > y1 and param2_i[i] < y3:
                param1.append(param1_i[i])
                param2.append(param2_i[i])
        param1 = np.array(param1)
        param2 = np.array(param2)

        xmin = param1.min()
        xmax = param1.max()
        ymin = param2.min()
        ymax = param2.max()
        xprime, yprime = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
        positions = np.vstack([xprime.ravel(), yprime.ravel()])
        values = np.vstack([param1, param2])
        kernel = stats.gaussian_kde(values)
        z = np.reshape(kernel(positions).T, xprime.shape)
        fig = plt.figure(figsize=(18,18))
        ax = fig.add_subplot(111)
        ax.imshow(np.rot90(z), cmap=plt.cm.gist_earth_r, extent=[xmin, xmax, ymin, ymax])
        ax.plot(param1, param2, 'k.', markersize=1)
        ax.set_xlim([xmin, xmax])
        ax.set_ylim([ymin, ymax])
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.title("Gaussian Kernel Density Estimation for " + ylabel + " vs " + xlabel)
        plt.savefig(os.path.join(self.directory, fname + '.png'))
        plt.show()
        r = np.corrcoef(param1, param2)[0][1]
        plt.figure(2)
        if diff:
            categories = self.searchhighestresids(75, colormap = True)
            categories = np.array(categories)
            colormap = np.array(['b', 'orange'])
            plt.scatter(param1_i, param2_i, s=0.75, c = colormap[categories])
            plt.xlabel(xlabel)
            plt.ylabel(ylabel)
            plt.title(ylabel + " vs " + xlabel + " Raw Data")
            plt.savefig(os.path.join(self.directory, fname2 + '.png'))
            plt.show()
        else:
            plt.scatter(param1_i, param2_i, s=0.75)
            plt.xlabel(xlabel)
            plt.ylabel(ylabel)
            plt.title(ylabel + " vs " + xlabel + " Raw Data")
            plt.savefig(os.path.join(self.directory, fname2 + '.png'))
            plt.show()
        return r, r**2

    def plot_3D(self, type , fname):
        x = self.data['XWIN']
        y = self.data['YWIN']
        oceta = self.data['O-Ceta']
        oxi = self.data['O-Cxi']
        resids = [(((oceta[i] ** 2) + (oxi[i] ** 2)) ** .5) for i in range(len(oceta))]
        if type == 'trisurf':
            z = resids
            fig = plt.figure()
            ax = fig.gca(projection='3d')
            ax.set_xlabel('XWIN')
            ax.set_ylabel('YWIN')
            ax.set_zlabel('magnitude of residual')
            ax.plot_trisurf(x, y, z, cmap='magma')
            plt.savefig(os.path.join(self.directory, fname + '.png'))
            plt.show()

        elif type == 'quiver':
            x = np.array(x)
            y = np.array(y)
            over = self.searchhighestresids(90)
            #x = [self.data['XWIN'][i] for i in over]
            #y = [self.data['YWIN'][i] for i in over]
            z = np.array(resids)
            u = np.array(self.data['O-Ceta'])
            v = np.array(self.data['O-Cxi'])
            w = np.array([0 for i in range(len(v))])
            p3.clear()
            quiver = p3.quiver(x,y,z,u,v,w, size = 3)
            #p3.savefig(os.path.join(self.directory, fname + '.png'))
            p3.show()

    def searchhighestresids(self, percentile, colormap = False): #percentile is a number between 0 and 100
        if not colormap:
            oceta = self.data['O-Ceta']
            oxi = self.data['O-Cxi']
            resids = [(((oceta[i] ** 2) + (oxi[i] ** 2)) ** .5) for i in range(len(oceta))]
            resid_thresh = np.percentile(resids, percentile)
            selected_indices = []
            for i in range(len(resids)):
                if resids[i] > resid_thresh:
                    selected_indices.append(i)
            return selected_indices
        else:
            oceta = self.data['O-Ceta']
            oxi = self.data['O-Cxi']
            resids = [(((oceta[i] ** 2) + (oxi[i] ** 2)) ** .5) for i in range(len(oceta))]
            resid_thresh = np.percentile(resids, percentile)
            categories = []
            for i in range(len(resids)):
                if resids[i] > resid_thresh:
                    categories.append(1)
                else:
                    categories.append(0)
            return categories

    def getresidpercentile(self, percentile, fwhm = False):
        if not fwhm:
            oceta = self.data['O-Ceta']
            oxi = self.data['O-Cxi']
            resids = [(((oceta[i] ** 2) + (oxi[i] ** 2)) ** .5) for i in range(len(oceta))]
            resid_thresh = np.percentile(resids, percentile)
            mindelta = 1000
            winner = 0
            for i in range(len(resids)):
                if np.abs(resid_thresh - resids[i]) < mindelta:
                    mindelta = np.abs(resid_thresh - resids[i])
                    winner = i
            return winner
        else: #FWHM/MagAuto
            resids = [self.data['FWHM'][i]/self.data['MagAuto'][i] for i in range(len(self.data['FWHM']))]
            resid_thresh = np.percentile(resids, percentile)
            mindelta = 1000
            winner = 0
            for i in range(len(resids)):
                if np.abs(resid_thresh - resids[i]) < mindelta:
                    mindelta = np.abs(resid_thresh - resids[i])
                    winner = i
            return winner

    def get_mean_subplot_value(self, star):
        image_data = fits.getdata(self.fits_file)
        bound = 20
        x = self.data['XIM'][star]
        y = self.data['YIM'][star]
        if y-bound < 0:
            ylow = 0
        else:
            ylow = int(y-bound)
        if x-bound < 0:
            xlow = 0
        else:
            xlow = int(x-bound)
        if y+bound > len(image_data):
            yup = len(image_data)
        else:
            yup = int(y+bound)
        if x+bound > len(image_data[0]):
            xup = len(image_data[0])
        else:
            xup = int(x+bound)
        image_data = image_data[ylow:yup, xlow:xup]
        pixelvals = np.array(image_data.flatten())
        meanval = np.mean(pixelvals)
        return meanval

    def get_stardata(self, star):
        image_data = fits.getdata(self.fits_file)
        bound = 20
        x = self.data['XIM'][star]
        y = self.data['YIM'][star]
        # print(image_data)
        # image_data = image_data[int(y)-20: int(y)+20][int(x)-20: int(x)+20]
        ylow, yup, xlow, xup = 0, 0, 0, 0
        if y - bound < 0:
            ylow = 0
        else:
            ylow = int(y - bound)
        if x - bound < 0:
            xlow = 0
        else:
            xlow = int(x - bound)
        if y + bound > len(image_data):
            yup = len(image_data)
        else:
            yup = int(y + bound)
        if x + bound > len(image_data[0]):
            xup = len(image_data[0])
        else:
            xup = int(x + bound)
        image_data = image_data[ylow:yup, xlow:xup]
        centroid = (self.data['XIM'][star]-xlow, self.data['YIM'][star]-ylow)
        centroid2 = (self.data['XWIN'][star]-xlow, self.data['YWIN'][star]-ylow)
        annulus = self.data['FWHM'][star]

        return image_data, annulus, centroid, centroid2

    def mpv_excluding_annulus(self, star):
        image_data, annulus, centroid, centroid2 = self.get_stardata(star)
        #exclude pixel values within annulus
        nap = 0
        pixcount = 0
        for i in range(len(image_data)):
            for j in range(len(image_data[0])):
                placeholder = [[j - 0.25, i - 0.25],
                               [j - 0.25, i + 0.25],
                               [j + 0.25, i - 0.25],
                               [j + 0.25, i + 0.25]]
                for k in placeholder:
                    if annulus ** 2 < (k[0] - centroid[0]) ** 2 + (k[1] - centroid[1]) ** 2:
                        pixcount += image_data[i][j] / 4
                        nap += 0.25
        mpv = pixcount/nap
        return mpv

    def starstatgraph(self, star): #star is an index in data
        lighthist_bins = 150
        residhist_bins = 150
        maghist_bins = 150
        image_data, annulus, centroid, centroid2 = self.get_stardata(star)

        #let the plots begin
        plt.figure(1)
        fig, ax = plt.subplots()
        plt.imshow(image_data, aspect="auto")
        plt.title("Star at Index " + str(star))
        circle = plt.Circle(centroid, annulus, fill = False, color = 'red')
        ax.add_artist(circle)
        plt.colorbar()
        plt.scatter(centroid[0],centroid[1], color='red', s = 1)
        #plt.scatter(centroid2[0], centroid2[1], color = 'orange', s = 0.5)
        plt.savefig(os.path.join(self.directory, str(star) + 'a.png'))
        plt.show()

        plt.figure(2)
        if not os.path.exists(os.path.join(self.directory, self.number + '_mpvs.csv')):
            mpv = pd.Series([self.mpv_excluding_annulus(star) for star in range(len(self.data['Eta']))])
            mpv.to_csv(os.path.join(self.directory, self.number + '_mpvs.csv'))

        mpv = np.array(pd.read_csv(os.path.join(self.directory, self.number + '_mpvs.csv'), header=None)[1])
        n, bins, patches = plt.hist(mpv,lighthist_bins)
        patchindex = 0
        for i in range(1, len(bins) - 1):
            if mpv[star] > bins[i - 1] and mpv[star] < bins[i + 1]:
                patchindex = i - 1
        patches[patchindex].set_fc('r')
        plt.title("Subplot Pixel Value Distribution")
        z = (mpv[star] - np.mean(mpv)) / np.std(mpv)
        zscore = "z = " + str(z)
        px = np.max(mpv) - 0.25 * (np.max(mpv) - np.min(mpv))
        py = np.max(n) - 0.25 * (np.max(n) - np.min(n))
        plt.text(px, py, zscore, ha='center', va='center')
        plt.savefig(os.path.join(self.directory, str(star) + 'b.png'))
        plt.show()

        plt.figure(3)
        oceta = self.data['O-Ceta']
        oxi = self.data['O-Cxi']
        resids = [(((oceta[i]**2)+(oxi[i]**2))**.5) for i in range(len(oceta))]
        n, bins, patches = plt.hist(resids, bins=residhist_bins, align='left', color='g')
        patchindex = 0
        for i in range(1,len(bins)-1):
            if resids[star] > bins[i-1] and resids[star] < bins[i+1]:
                patchindex = i-1
        patches[patchindex].set_fc('r')
        plt.title("Distribution of Residual Magnitudes")
        z = (resids[star] - np.mean(resids))/np.std(resids)
        zscore = "z = " + str(z)
        px = np.max(resids) - 0.25 * (np.max(resids) - np.min(resids))
        py = np.max(n) - 0.25 * (np.max(n) - np.min(n))
        plt.text(px, py, zscore, ha='center', va='center')
        plt.savefig(os.path.join(self.directory, str(star) + 'c.png'))
        plt.show()

        plt.figure(4)
        mag = self.data['MagAuto']
        n, bins, patches = plt.hist(mag, bins=maghist_bins, align='left', color='y')
        patchindex = 0
        for i in range(1, len(bins) - 1):
            if mag[star] > bins[i - 1] and mag[star] < bins[i + 1]:
                patchindex = i - 1
        patches[patchindex].set_fc('r')
        plt.title("Distribution of Stellar Magnitudes")
        z = (mag[star] - np.mean(mag)) / np.std(mag)
        zscore = "z = " + str(z)
        px = np.min(mag) + 0.25*(np.max(mag) - np.min(mag))
        py = np.max(n) - 0.25*(np.max(n) - np.min(n))
        plt.text(px, py, zscore, ha='center', va='center')
        plt.savefig(os.path.join(self.directory, str(star) + 'd.png'))
        plt.show()




#starfield = StarField('20190823', '0357', data, fits_file)
#starfield.write_outfile()