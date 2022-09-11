#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  8 11:30:18 2022

@author: sergio
"""
# (and the mean values)
import numpy as np
from tqdm import tqdm 
from astropy.io import fits
from healpy.fitsfunc import read_map 
import os
from unit_converter import transform, MJy2K



# def save(data_origin,frecs, direc='./data/', data_file='data.fits', data_file_mean='mean.fits'):

    
# data_mean: str. Default is 'mean.fits'
#         Name of the file where the data means are going to be stored
def save(data_origin, frecs, direc='../../output/simulation_1/raw_data/', data_file='data_raw.fits', sum_file='suma_raw.fits'):
    """                   
    Given several .fits data files, mean for each pixel is computed.
    It saves the data in a .fits file 
    
    Parameters
    ----------
    data_origin: List
        paths where all the data files are stored
    
    direc: str. Default is ./output/raw/
        Folder where the data is going to be stored
    
    data_file: str. Default is 'data_raw.fits'
        Name of the file where all the data is going to be stored
    
    
    
    """  
    # create folder if needed
    if not os.path.exists(direc):
        os.mkdir(direc)
    f1 = str(direc)+str(data_file)
    f2 = str(direc)+str(sum_file)
    
    if os.path.exists(f1):
        if input('Files already found. Overwrite? [y/n]: ') != 'y':
            print('aborted')
            return
        if os.path.exists(f1):
            os.remove(f1)
        if os.path.exists(f2):
            os.remove(f2)
    #Check if files already exists, if so, ask for permission to overwrite
    
    
    #First, read the files
    
    
    #Empty lists to add the data and the units
    data = []
    #Go down the raw data files
    for j, i in enumerate(tqdm(data_origin)):
        data.append(fits.getdata(i))
        #Check units
        # hdr = fits.open(i)[1].header
        # units = hdr['TUNIT1']
        # angle = hdr['BEAMSIZE']
        # if units == 'MJy/sr':
        #     aux = read_map(i, dtype=np.float64)
        #     data.append(MJy2K(aux, float(frecs[j])))
        #     # angle = hdr['BEAMSIZE']
        #     # factor = factor_MJy(frecs[j], angle)
        #     # data.append(read_map(i, dtype=np.float64)*factor)
        # else:
        #     factor = transform[units]
        #     #Read the file, transform to K and add to data list
        #     data.append(read_map(i, dtype=np.float64)*factor)
        # aux = fits.open(i)[1].data.field(0)
        # data.append(np.ndarray.flatten(aux.astype(np.float64))*factor)
    #Convert to array for treatment
    data = np.array(data).T
    
    # Compute sum
    suma = np.sum(data, axis=1)
    
    #Save data files
    h1 = fits.HDUList([fits.PrimaryHDU(data)])
    h2 = fits.HDUList([fits.PrimaryHDU(suma)])
    h1.writeto(f1)
    h2.writeto(f2)
    
    
            
    # return data
        
    
    
    
    