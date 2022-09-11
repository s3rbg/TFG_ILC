#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 24 09:48:37 2022

@author: sergio
"""

from healpy.fitsfunc import read_map
import healpy as hp
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import tri
from scipy.optimize import curve_fit as cf
from astropy.io import fits
from uncertainties import ufloat
from scipy.interpolate import griddata
def f(x, a, b):
    return a*x + b

def fit(f, x, y):
    popt, pcov = cf(f, x, y)
    return popt, np.sqrt(np.diag(pcov))

def fit2(x, y):
    popt, pcov = np.polyfit(x, y, deg=1, cov=True)
    return popt, np.sqrt(np.diag(pcov))
    
def tt_plot(cmb, result, mini, maxi, title, mask=np.array([0])):
    if len(mask)==1:
        popt, pcov = fit(f, cmb, result)
        plt.plot(cmb, result, 'b.')
    else:
        popt, pcov = fit(f, cmb[mask], result[mask])
        plt.plot(cmb[mask], result[mask], 'b.')
    x = np.linspace(mini, maxi, 1000)
    y = f(x, popt[0], popt[1])  
    plt.plot(x, y, 'k-')
    plt.axis('equal')
    plt.xlabel(r'$T(CMB)$/$\mu$K')
    plt.ylabel(r'$T(result)$/$\mu$K')
    plt.xlim((mini, maxi))
    plt.ylim((mini, maxi))
    plt.title(title)
    plt.grid()
    plt.show()
    a = ufloat(popt[0], pcov[0]) 
    b = ufloat(popt[1], pcov[1])
    return a, b

# file = '../../CMB_map.fits'
cmb = fits.getdata('../../CMB/mean.fits')*1e6
hp.mollview(cmb, norm='hist', title='CMB', unit='uK')
# cmb = read_map('../../Components/CMB/map_100.fits')*1e6
# hp.mollview(cmb, norm='hist', title='CMB', unit='uK')
# cmb = hp.pixelfunc.ud_grade(cmb, 1024, order_in='RING', order_out='RING')
# cmb = fits.getdata(file)*1e6
# hp.mollview(cmb, norm='hist', title='CMB')
# plt.show()
suma_reg = fits.getdata('../../output/simulation_1/output_data/suma_fin.fits')*1e6
suma_full = fits.getdata('../../output/simulation_2/output_data/suma_fin.fits')*1e6

mask = fits.getdata('../../mask.fits').astype(np.bool_)
# mask = fits.getdata('../../Masks/mask40.fits').astype(np.bool_)
hp.mollview(suma_full, title='Full Sky', norm='hist', unit='uK')
hp.mollview(suma_reg, title='Regions', norm='hist', unit='uK')
suma_full_mask = hp.ma(suma_full)
suma_full_mask.mask = np.logical_not(mask)
hp.mollview(suma_full_mask, title='Full Sky', norm='hist', unit='uK')

suma_reg_mask = hp.ma(suma_reg)
suma_reg_mask.mask = np.logical_not(mask)
hp.mollview(suma_reg_mask, title='Regions', norm='hist', unit='uK')

cmb_mask = hp.ma(cmb)
cmb_mask.mask = np.logical_not(mask) 
hp.mollview(cmb_mask, title='CMB', norm='hist')
plt.show()



a_reg, b_reg = tt_plot(cmb, suma_reg, -600, 600, 'Regions', mask)
a_full, b_full = tt_plot(cmb, suma_full, -600, 600, 'Full Sky', mask)

print('Regions y=ax+b:')
print('a =', a_reg)
print('b =', b_reg)

print('\nFull sky y=ax+b:')
print('a =', a_full)
print('b =', b_full)

