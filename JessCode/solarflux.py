#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 21 10:22:21 2021

@author: jesscmiel
"""

import numpy as np

# Solar flux function
# Equation 1 in Gough et. al

def sol_flux(t, star, planet):

    """ Calculates solar flux as a function of time (using mean distance of planet from Sun)

    Inputs:
    t: time [s] (equation only valid for t<t_Sun)

    Outputs:
    Solar flux [W/m^2]

    """

    L = ((1 + (2/5)*(1-t/4.5e9))**(-1))*(star.luminosity) # [W]
    S = L/(4*np.pi*(planet.dist**2)) # [W/m^2]
    return S
