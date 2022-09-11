#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 20 10:39:45 2022

@author: sergio
"""

from astropy.io import fits
import numpy as np
import healpy as hp
from healpy.fitsfunc import read_map
from unit_converter import MJy2K


# file = 'febecop_ffp10_lensed_scl_cmb_100_mc_0000.fits'
# frec = '030'
# file = 'psm_output/output/observations/LFI/detector_F'+str(frec)+'/group2_map_detector_F'+str(frec)+'.fits'
def smooth(target, frec, FWHM, save, obser=True):
    """
    Merges the maps after the smoothing process

    Parameters
    ----------
    target : float
        target FWHM in deg
    FWHM : float
        Resolution of the input map, in arcmin
    frec : float
        observation frecuency
    save : str
        saving path
    obser : bool, optional
        If False, only CMB is done. Otherwhise, the full observation
        map


    """
    target *= np.pi/180
    FWHM *= np.pi/180/60
    
    if obser:
        components = ['CMB/', 'Synchrotron/', 'Dust/', 'Freefree/']
    else:
        components = ['CMB/']
    # components = ['CMB', 'Dust', 'Freefree', 'Synchrotron']
    # save = 'Components/suma.fits'
    # frec = [30, 44, 70, 100, 143, 217, 353, 545, 857]
    # FWHM = [32.2879982, 26.99699974, 13.21800041, 9.862, 7.303, 5.021, 4.944, 4.831, 4.638]*np.pi/60/180
    
    aux = 'Components/'
    suma = np.zeros(12*1024**2)
    for i in components:
        file = str(aux)+str(i)+'map_'+str(frec)+'.fits'
        hdr = fits.open(file)[1].header
        unit = hdr['TUNIT1']
        nside = hdr['NSIDE']
        order = hdr['ORDERING']
        mapa = read_map(file)
        if nside == 2048:
            mapa = hp.pixelfunc.ud_grade(mapa, 1024, order_in=order, order_out='RING')
        
        
        x = np.sqrt(target**2-FWHM**2)
        mapa = hp.smoothing(mapa, x)
        if unit == 'MJy/sr':
            mapa = MJy2K(mapa, frec)
        suma+=mapa
    h1 = fits.HDUList([fits.PrimaryHDU(suma)])
    h1.writeto(save)

# def do_mean():

frecs = [30, 44, 70, 100, 143, 217, 353, 545, 857]
FWHM = [32.2879982, 26.99699974, 13.21800041, 9.862, 7.303, 5.021, 4.944, 4.831, 4.638]


def save_CMB():
    """
    Smoothes the CMB maps. 
    The files must be named map_{frec}, and stored in ./Components/CMB/.
    Make sure to have the destiny folder created as ./CMB/

    """
    #Define the saving folder (CMB/) and the beggining of the file
    folder = 'CMB/map_'
    target = 1 #Smooth up to 1 deg
    for i in range(9):
        file = str(folder)+str(frecs[i])+'.fits'
        smooth(target, frecs[i], FWHM[i], file, False)

def make_simulation():
    """
    Makes the smoothed simulation map.
    The files must be named map_{frec}, and stored in ./Components/{component name}/.
    {Component name} = CMB/, Dust/, Freefree/, Synchrotron/
    Make sure to have the destiny folder created as ./observations/
    """
    folder = 'observations/map_'
    target = 1 #Smooth up to 1 deg
    for i in range(9):
        file = str(folder)+str(frecs[i])+'.fits'
        smooth(target, frecs[i], FWHM[i], file, False)

def do_mean():
    """
    After having smoothed the CMB maps, compute its mean.
    The files must be stored in ./CMB, where the file "mean.fits"
    will be stored
    """
    folder = 'CMB/map_'
    suma = np.zeros(12*1024**2)
    for i in frecs:
        file = str(folder)+str(i)+'.fits'
        suma += fits.getdata(file)
    suma/=9
    h1 = fits.HDUList([fits.PrimaryHDU(suma)])
    h1.writeto('CMB/mean.fits')
    

    