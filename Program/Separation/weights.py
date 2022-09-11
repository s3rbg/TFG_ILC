#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  9 11:01:07 2022

@author: sergio
"""
import numpy as np
from astropy.io import fits
from tqdm import tqdm
import os

def get_weights(C):
    #Convert to matrix and invert
    C_m = np.asmatrix(C)
    C_inv = np.linalg.inv(C_m)
    
    #Denominator term
    den = np.sum(C_inv)
    
    #Column matrix of 1's, to avoid loops
    aux = np.matrix(np.full(len(C), 1)).T
    
    #Compute actual weight
    w = np.array(C_inv*aux/den).T[0]
    return w
        
def save_weights(C, folder='../../output/simulation_1/operators/', file = 'weight.fits'):
    #Create folder if necesary (just in case, it should have been created by now)
    if not os.path.exists(folder):
        os.mkdir(folder)
    #Get file path
    f1 = str(folder)+str(file)
    
    #Check if files already exist
    if os.path.exists(f1):
        if input('Weight files already found. Overwrite? [y/n]: ') != 'y':
            print('aborted')
            return
        os.remove(f1)
        
    #make C a 3D array (for coding simplicity)
    if len(np.shape(C))==2:
        C = np.array([C])
        
    #Initialize w    
    w = []
    
    #Compute weights for each covariance matrix ()
    for i in C:
        w.append(get_weights(i))
    #Save results
    h1 = fits.HDUList([fits.PrimaryHDU(w)])
    h1.writeto(f1)
    return np.array(w)

