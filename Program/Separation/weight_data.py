# -*- coding: utf-8 -*-
"""
Created on Thu Aug 11 13:19:59 2022

@author: Sergio
"""
import numpy as np
from tqdm import tqdm
from astropy.io import fits
from save_data import save
from evaluate import evaluate
import os


def weight_data(data, w, template=None, folder='../../output/simulation_1/output_data/', data_file='data_fin.fits', sum_file='suma_fin.fits'):
    """
    Given the weights over the data, computed using a template,
    it applies those weights and gets the final map

    Parameters
    ----------
    data : array-like
        Input data, in column format
    w : array-like
        weights of all the maps (each row is a region)
    template : array-like
        template used to compute the weights, with progressively 
        increasing int index, starting at 0

    """
    #Create folder if necesary
    if not os.path.exists(folder):
        os.mkdir(folder)
    # Get paths
    f1 = str(folder)+str(data_file)
    f2 = str(folder)+str(sum_file)
    # Check if existing files
    if os.path.exists(f1) or os.path.exists(f2):
        if(input('Final maps files found, overwrite? [y/n]: ')) != 'y':
            print('Aborted')
            return
        else:
            if os.path.exists(f1):
                os.remove(f1)
            if os.path.exists(f2):
                os.remove(f2)
    #If the full map at once, make w 1D->2D
    
    #Get regions values
    template, val = evaluate(template, len(data))
    
    #Evaluate in the regions 
    for i,j in enumerate(val):
        index = np.where(template==j)[0]
        datai = data[index]
        wi = w[i]
        data[index] *= w[i]
    
    # Sum to get the final maps
    suma = np.sum(data, axis=1)
    
    #Save final data
    h1 = fits.HDUList([fits.PrimaryHDU(data)])
    h1.writeto(f1)
    h1 = fits.HDUList([fits.PrimaryHDU(suma)])
    h1.writeto(f2)
# def weight_data(data, w, data_raw):
#     #Case 1 - full map, straight forward computation
#     if len(w.shape)==1:
#         data_fin = np.copy(data)
#         for i in range(len(w)):
#             data_fin[:,i] = data[:,i] * w[i]
#             #data ready for storage
#     #Case 2 - work with slices
#     else:
#         #   Initialize final array as a line of 0s
# ####Note: Remember to erase that line at the end
        
#         data_fin = np.array([[0]*len(w[0])])
        
       
        
#         for i in tqdm(range(len(w))):
#             data_i = data[i]
#             data_aux = np.copy(data_i)
#             w_i = w[i]
#             for j in range(len(w_i)):
#                 data_aux[:,j] = data_i[:,j]*w_i[j]
#             #slice of data weighted, store in data final
#             data_fin = np.concatenate((data_fin, data_aux))
#         data_fin = np.delete(data_fin, 0, 0)
#         data_fin=np.delete(data_fin, slice(len(data_raw), len(data_fin)), 0)
#         map_fin = []
#         for i in tqdm(data_fin):
#             map_fin.append(np.sum(i))
#         h1 = fits.HDUList([fits.PrimaryHDU(data_fin)])
#         h2 = fits.HDUList([fits.PrimaryHDU(map_fin)])
#         h1.writeto('output/final/data_w.fits')
#         h2.writeto('output/final/suma.fits')