#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  8 12:54:37 2022

@author: sergio
"""
import numpy as np
#Dictionary to convert data to kelvin
transform = {'K_CMB': 1, 'uK_CMB': 1e-6}

T0 = 2.7255#K, mean temperature of the CMB
c = 3e8 #m/s, speed of light
k_B = 1.380649e-23 #J/K, Boltzmann constant
h = 6.62607015e-34 #J s
factor = 1.472e-3 #MJy/sr, factor that appears in Plancks law when 
# frequency is given in GHz (2h/c²)*(1GHz)³


def MJy2K(data, frec):
    """
    
    #Conversion factor from MJy/sr to K
    Parameters
    ----------
    frec : float, string
        frecuency, in GHz

    Returns
    -------
    Factor to convert to K

    """
    frec = float(frec)*1e9
    t2 = (c**2*k_B*T0**2)/(2*h**2*frec**4)
    exp = np.exp((h*frec)/(k_B*T0))
    num = (exp-1)**2
    return data*num/exp*t2*1e-20
    
    
    
