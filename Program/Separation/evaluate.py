#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 20 09:32:05 2022

@author: sergio
"""
import numpy as np
def evaluate(template=None, n_pix=12*1024**2):
    """
    Gets the values used to separate the regions for the future evaluation
    of the map

    Parameters
    ----------
    template : array-like
        Separation templation, with the same number of pixels as the 
        maps and not repeated indexes in separated regions. 
        If not given, taken all map by default
        
    n_pix: int
        number of pixels in the map, needed if no template is 
    given
    
    Returns
    ------- 
    template: array-like 
        if a template is given, this is a copy. Otherwise, one is created
        with only one region covering the full map
    Values: array-like

    """
    #Get the list of different values used in the template
    if isinstance(template, np.ndarray):
        val = list(set(template))
    #If no template, save as 0s array (all map has the 0 index, so no regions are separated)
    elif isinstance(template, type(None)):
        template = np.zeros(n_pix)
        val = [0]
    return template, val