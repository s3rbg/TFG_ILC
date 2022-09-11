#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 19 09:36:30 2022

@author: sergio
"""
import numpy as np
from astropy.io import fits
from healpy.pixelfunc import pix2ang
import healpy as hp
import os
def make_template(lat, nside=1024, ipix=12*1024**2, nest = False):
    """
    Divides the sky map in regions, creating a template

    Parameters
    ----------
    lat : float or array-like
        latitude to separate at, in deg
    nside : int or array-like
      The healpix nside parameter, must be a power of 2, less than 2**30
    ipix : int or array-like
      Pixel indices
    nest : bool, optional
      if False, assume RING pixel ordering, otherwise, NESTED pixel ordering

    """
    if os.path.exists('output/aux/template.fits'):
        if(input('File found, overwrite? [y/n]: ')) != 'y':
            print('Aborted')
            return
        else:
            os.remove('output/aux/template.fits')
    lat = np.radians(lat)
    index = np.arange(0, ipix)#Array of indexes
    temp = np.zeros(ipix) #Template, initialized as 0
    lati, long = pix2ang(nside, index) #ang
    lati -= np.pi/2 #Change reference to 0 at ecuator for easier computation
    
    #Convert input latitude to array with negative and positive values (regions delimiters)
    if type(lat) == list:
        lat = np.array(lat)
    if type(lat) == int:
        lat = float(lat)
        
    if type(lat)==np.float64:
        lat = np.array([-lat, lat])
    else:
        lat = lat.astype(float)
        aux = np.flip(-lat)
        lat = np.concatenate((aux, lat))
        
    #Start actual delimitation
    n_it = len(lat)+1
    for i in range(n_it):
        #Borders
        #Upper part
        if i == 0:
            temp[np.where(lati<lat[i])[0]] = i+1
        #Lower part
        elif i==n_it-1:
            temp[np.where(lati>lat[i-1])] = 1
        #Middle
        else:
            #Consider all the regions but the ecuatorial
            if np.sign(lat[i-1]) == np.sign(lat[i]):
                #Criteria is that borderline goes in the region closer to the ecuator
                if lat[i]<0:
                    temp[np.where((lati>=lat[i-1])&(lati<lat[i]))] = i+1
                else:
                    temp[np.where((lati>lat[i-1])&(lati<=lat[i]))] = i
                # aux = np.abs([lat[i-1], lat[i]])
                # temp[np.where((np.abs(lati)<=max(aux))&(np.abs(lati)>min(aux)))] = i+1
        
    h1 = fits.HDUList([fits.PrimaryHDU(temp)])
    h1.writeto('templatev8.fits')
make_template([25])