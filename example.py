from utilities import *

#download files from pipeline:
#fetchlmi('20190901a', '29P', 'Rb')
#read in example starfield:
data, image_data = setuplmi(rpath='29p-data/20190901a/R/20190901a.R.0147',
                            rbpath = '29p-data/20190901a/Rb/20190901a.Rb.0147')

#analyze example starfield:
"""data, image_data = setuplmi(rpath='29p-data/20190901a/R/20190901a.R.0147',
                            rbpath = '29p-data/20190901a/Rb/20190901a.Rb.0147',
                            make_output=True)"""
#graphically search for star data:
findstar(rpath='29p-data/20190901a/R/20190901a.R.0147',
        rbpath = '29p-data/20190901a/Rb/20190901a.Rb.0147',
        destdir = '29p-data/')
