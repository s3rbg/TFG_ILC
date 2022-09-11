#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 25 11:37:39 2022

@author: sergio
"""
import healpy as hp
from healpy.fitsfunc import read_map
from astropy.io import fits
import numpy as np


def reminder(cmb, mapa, title):
    rem_map = mapa-cmb
    hp.mollview(rem_map, title=title, norm='hist', unit='uK')
    # hp.mollview(rem_map, title=title, unit='uK')
    print('Mean '+str(title)+' = '+str(np.mean(rem_map)))
    print('StDev '+str(title)+' = '+str(np.std(rem_map)))


cmb = fits.getdata('../../CMB/mean.fits')*1e6
# cmb = read_map('Mapa_CMB.fits', field=0, dtype=np.float64, nest=True)
# name = 'PLA_simulation'
# cmb = hp.pixelfunc.ud_grade(cmb, 1024, order_in='NESTED', order_out='RING')*1e6

map_full = fits.getdata('../../output/simulation_2/output_data/suma_fin.fits')*1e6
map_region = fits.getdata('../../output/simulation_1/output_data/suma_fin.fits')*1e6

mask = fits.getdata('../../mask.fits').astype(np.bool_)

cmb_mask = hp.ma(cmb)
cmb_mask.mask = np.logical_not(mask)

map_full_mask = hp.ma(map_full)
map_full_mask.mask = np.logical_not(mask)

map_region_mask = hp.ma(map_region)
map_region_mask.mask = np.logical_not(mask)
print('Intrinsec map parameters\nNo mask\n\nCMB')
print('Mean:',np.mean(cmb))
print('STD:',np.std(cmb))
print('\nFull Sky')
print('Mean',np.mean(map_full))
print('STD',np.std(map_full))
print('\nRegions')
print('Mean',np.mean(map_region))
print('STD',np.std(map_region))
print('\n\nWith mask\n\nCMB')
print('Mean',np.mean(cmb_mask))
print('STD',np.std(cmb_mask))
print('\nFull Sky')
print('Mean',np.mean(map_full_mask))
print('STD',np.std(map_full_mask))
print('\nRegions')
print('Mean',np.mean(map_region_mask))
print('STD',np.std(map_region_mask))


print('Results of the reminder plot\n\nWithout mask:')
reminder(cmb, map_full, 'Full Sky')
reminder(cmb, map_region, 'Regions')
print('\nWith mask:')
reminder(cmb_mask, map_full_mask, 'Full Sky')
reminder(cmb_mask, map_region_mask, 'Regions')
