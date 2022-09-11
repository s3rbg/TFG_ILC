# -*- coding: utf-8 , -*-
"""
Created on Thu Aug 11 10:35:51 2022

@author: Sergio
"""

import numpy as np
from astropy.io import fits
from save_data import save
import numpy.lib.recfunctions as rfn
import pandas as pd
from covariance_regions import covar
from weights import save_weights
from weight_data import weight_data
# def main():
#     path = ['group2_map_detector_F070.fits']
#     save(path, None)
def merge():
    master_path = '../../Observations/map'
    
    frecs = [30, 44, 70, 100, 143, 217, 353, 545, 857]
    # frecs_high = ['100', '143', '217', '353']
    paths = []
    for i  in frecs:
        paths.append(str(master_path)+str(i)+'.fits')
    save(paths, frecs, '../../output/raw_data/')    
# merge()

def do_cov():
    data = fits.getdata('../../output/raw_data/data_raw.fits')
    # mean = fits.getdata('output/raw/mean_raw.fits')
    template = fits.getdata('../../templatev6.fits')
    print('Doing cov')
    covar(data, template, '../../output/simulation_1/operators/')

def do_weight():
    C = fits.getdata('../../output/simulation_1/operators/cov.fits')
    print('Doing weights')
    save_weights(C, '../../output/simulation_1/operators/')
    
def comp_weight():
    template = fits.getdata('../../templatev6.fits')
    w = fits.getdata('../../output/simulation_1/operators/weight.fits')
    data = fits.getdata('../../output/raw_data/data_raw.fits')
    print('\nWeighting map')
    weight_data(data, w, template, '../../output/simulation_1/output_data/')
    
# merge()
do_cov()
do_weight()
comp_weight()
# 
