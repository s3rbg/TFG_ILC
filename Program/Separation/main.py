#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  9 18:01:54 2022

@author: sergio
"""

from save_data import save
from astropy.io import fits
from covariance_regions import covar
from weights import save_weights
from weight_data import weight_data
def prepare_files(fold_origin='../../Observations/', fold_destiny='../../output/raw_data'):
    """
    Given the folder where the observationrs are, prepare them for the rest
    of the program. The files must be named 'map{frecuency}.fits', where
    {frequency} is the frequency of the observation, in GHz

    Parameters
    ----------
    fold_origin : str
        folder wrehe the observaitons are
    fold_destiny : str
        folder where the raw data, prepared to be applied the code,
        is stored

    """
    #List of frecuencies
    frecs = [30, 44, 70, 100, 143, 217, 353, 545]
    #Initialize list of file paths
    paths = []
    #Fill the list
    for i  in frecs:
        paths.append(str(fold_origin)+'map_'+str(i)+'.fits')
    save(paths, frecs, fold_destiny)
    
def separate(data_file='../../output/raw_data/data_raw.fits', temp_file='None', 
             destiny = '../../output/simulation_1/' ):
    """
    Performs the separation process

    Parameters
    ----------
    data_file : str
        Path of the 'data_raw.fits' file resulted from prepare_files()
    temp_file : str
        Path of the template of separation used. If 'None', full sky will be taken
    destiny : str
        destiny folder
    """
    #Open files
    data = fits.getdata(data_file)
    #Open template, if given
    if temp_file != 'None':
        template = fits.getdata(temp_file)
    else:
        template = None
    
    #Operators folder, where C and w will be stored
    op_folder = str(destiny)+'operators/'
    # Compute covariance
    C = covar(data, template, op_folder)
    print('Covariance done')
    # Compute weights
    w = save_weights(C, op_folder)
    print('Weights done')
    
    #result folder
    fin_folder = str(destiny)+'output_data/'
    weight_data(data, w, template, fin_folder)

# prepare_files('../../observations/' ,'../../output/raw_data/')
separate('../../output/raw_data/data_raw.fits' , '../../templatev2.fits', destiny = '../../output/simulation_4/')
# separate('../../output/raw_data/data_raw.fits' , destiny = '../../output/simulation_2/')   