#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 21 10:29:32 2021

@author: jesscmiel
"""

import numpy as np
import meltfraction

# Function to solve for roots = P_x
def func(x, N_x_tot, planet, star, gas):
    """Define idealized expression with x representing pressure so that pressure can be solved for (using Newton's method).
    
    Inputs: 
        
    x: pressure [bar]
    N_x_tot: total number of moles in system
    t: time [yr]
    planet
    star: star which planet orbits
    gas: species of gas
    
    Outputs:
    f: idealized expression to be evaluated
    """
    
    #ppmw to total mass fraction
    ppmw_convert = 1e-6
    
    # #mass fraction of magma ocean 
    # f_mo = 0.0001
    
    
    
    # Defining constants to idealize expression
    B = (4*np.pi*planet.radius**2)/(planet.g*gas.mu)
    C = ((planet.mass*(1-planet.CMF))*(gas.alpha*ppmw_convert))/gas.mu
    
    f = B*x + C*meltfraction.psi(x, planet, star, gas)*x**(1/gas.beta) - N_x_tot
    return f
