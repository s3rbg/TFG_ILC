# -*- coding: utf-8 -*-
"""
Created on Tue Aug 16 22:09:05 2022

@author: Sergio
"""

import numpy as np
from astropy.io import fits
from tqdm import tqdm
import os



def cov(maps):
    """
    Computes the covariance matrix
    of the given maps

    Parameters
    ----------
    maps : array-like
        raw data of the maps, in column format
    mean : array-like
        data with the mean value for each pixel

    Returns
    -------
    Covariance matrix

    """
    
    
    n_pix = len(maps)
    #Initialize matrix as an array
    maps = np.matrix(maps)
    
    return 1/n_pix*maps.transpose()*maps
    