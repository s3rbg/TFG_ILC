#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 19 19:45:39 2022

@author: sergio
"""
import numpy as np
from evaluate import evaluate
from astropy.io import fits
import healpy as hp 
from covariance2 import cov
import os

# def covar(maps, template=None):
def covar(maps, template=None, folder='../../output/simulation_1/operators/', file_name='cov.fits'):
    """
    

    Parameters
    ----------
    maps : array-like
        column format array with the data of all the maps
    template : array-like
        Separation templation, with the same number of pixels as the 
        maps and not repeated indexes in separated regions. 
        If not given, taken all map by default
    folder:str
        destiny folder
    """
    #Create folder if necesary
    if not os.path.exists(folder):
        os.mkdir(folder)
    #get full path
    file = str(folder)+str(file_name)
    #Check for existing files
    if os.path.exists(file):
        if(input('Covariance file found, overwrite? [y/n]: ')) != 'y':
            print('Aborted')
            return
        else:
            os.remove(file)
    template, val = evaluate(template, len(maps))
    C = []
    #Loop over all regions
    for i in val:
        #Compute covariance matrix in region "i"
        index = np.where(template==i)[0]
        C.append(cov(maps[index]))
    
    #Save in .fits file
    h1 = fits.HDUList([fits.PrimaryHDU(C)])
    h1.writeto(file) 
    
    return np.array(C)
    