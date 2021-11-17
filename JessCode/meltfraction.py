#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 21 10:26:41 2021

@author: jesscmiel
"""
from scipy import special
#import solarflux
import surfacetemp

# Mass of magma ocean as a function of surface temperature
def psi(P_s, planet, star, gas):
    """Calculate total silicate layer melt fraction Ψ.
    
    Inputs:
    P_s: Surface pressure [Pa] 
    planet: planet
    star: star that planet orbits
    gas: species of gas
    
    Outputs:
    melt_frac: melt fraction of the silicate layer [nondim]
    """
    melt_frac = .5*(special.erf((surfacetemp.surf_temp(P_s, star.t, planet, star, gas) - planet.T_t)/planet.dT)+1)
    return melt_frac
    