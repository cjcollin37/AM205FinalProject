#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 21 10:25:04 2021

@author: jesscmiel
"""
import solarflux

# Stefan-Boltzmann Constant
sigma = 5.670373e-8 # [W/(m^2 K^4)]

# Surface temperature function
def surf_temp(P_s, t, planet, star, gas):
    """Calculate surface temperature as a function of surface pressure over time.

    Inputs:
    P_s: surface pressure [Pa]
    t: time [yr]
    planet: planet
    star: star that planet orbits
    gas: species of gas

    Outputs:
    Surface temperature [K] 
    """
    # Energy absorbed by planet
    F = (1/4)*(1-planet.A)*solarflux.sol_flux(t, star, planet) #[W/m^2]
    # Equilibrium temperature 
    # Energy emitted (Stefan Boltzmann Law) = energy absorbed solved for equilibrium temp
    T_eq = (F / sigma)**(1/4) # [K]
    
    # Equation for surface temperature
    T_s = T_eq*(planet.P_em/P_s)**(-gas.R/gas.c_p) # K
    
    # Print formatting
    # print('Surface temperature: {} K'.format(T_s))
    return T_s

