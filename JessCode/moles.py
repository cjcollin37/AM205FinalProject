#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 21 10:27:52 2021

@author: jesscmiel
"""

import numpy as np
import meltfraction

#ppmw to total mass fraction
ppmw_convert = 1e-6

def partitioning(P_s, planet, star, gas):
    """Calculate moles of species present in magma ocean and atmosphere.
   
    Inputs:
    P_s: surface pressure [Pa]
    planet: planet
    star: star that planet orbits
    gas: species of gas

    
    Outputs:
    N_x_atm: total number of moles of species in atmosphere
    N_x_mo: total number of moles of species in magma ocean
    """
    # Calculate mass of species in atmosphere from surface pressure
    M_atm = ((P_s)*(4*np.pi*(planet.radius**2)))/(planet.g) #[kg]
    
    # Calculate moles of species in atmosphere
    N_atm = M_atm / gas.mu # [kg/mol]
    
    # Mass of magma ocean
    M_mo = planet.mass*(1-planet.CMF)*meltfraction.psi(P_s, planet, star, gas)

    # Mass of species in magma
    M_x_mo = M_mo * gas.alpha*ppmw_convert*((P_s)**(1/gas.beta))

    # Calculate moles of CO2 in MO
    N_x_mo = M_x_mo / gas.mu # [kg/mol]

    print(N_atm, ' moles in atmosphere')
    print(N_x_mo, ' moles in magma ocean')
    return N_atm, N_x_mo