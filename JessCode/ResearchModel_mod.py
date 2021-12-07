#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 10 14:37:23 2021

@author: jesscmiel
"""

#   To-do:
# ------------
# - Figure out how to format plot labels with f-string and latex subscripts

# - Validate the code! Make sure everything is working correctly and conserving moles of species
# - Change volume of magma ocean, see what partitioning is- check that you get constant for N_tot
#       - Plot melt fraction (psi) (define array for psi in the model) vs. moles: MO, atm., total
#       - Produce this plot for CO2, then for different species to see how it varies
# - Chemical dependence of solubility





# Modeling Magma Ocean and Atmospheric Interaction on Venus
# Import relevant libraries
import numpy as np
import matplotlib.pyplot as plt
import scipy.misc
import sys


# Import modules
import func
import func_mod
import moles


# Define classes
class Planet:
    """
    Defines the characteristics of a planet.
    Parameters:
    -----------
    mass: mass [kg]
    radius: radius of planet[m]
    g: Gravitational acceleration on planet [m/s^2]
    A: albedo [nondim]
    CMF: core mass fraction of planet [nondim]
    dist: distance from planet to Sun [m]
    T_t: [K]
    dT: [K]
    P_em: emitted pressure [Pa]
    """
    
    def __init__(self, mass, radius, g, A, CMF, dist, T_t, dT, P_em):
       self.mass = mass
       self.radius = radius
       self.g = g
       self.A = A
       self.CMF = CMF
       self.dist = dist
       self.T_t = T_t
       self.dT = dT
       self.P_em = P_em 
       
class Gas:
    """
    Defines the characteristics of a gas.
    
    Parameters:
    -----------
    R: Specific gas constant [J/(kgK)]
    c_p: Isobaric specific heat capacity [J/(kgK)]
    mu: molar mass of species [kg/mol]
    alpha: Partitioning constant for species [ppmw/Pa] from Lichtenberg 
        et. al 2021
    beta: Henry coefficient parameter [nondim] from Lichtenberg et. al 2021
    """

    def __init__(self, R, c_p, mu, alpha, beta):
        self.R = R
        self.c_p = c_p
        self.mu = mu
        self.alpha = alpha
        self.beta = beta

class Star:
    """
    Defines characteristics of a star.
    
    Parameters:
    -----------
    luminosity: solar luminosity [W]
    t: time [yr]
    """
    
    def __init__(self, luminosity, t):
        self.luminosity = luminosity
        self.t = t



# Newton's method function to solve for pressure
def newton(x0, N_x_tot, planet, star, gas, melt_fraction, tol=1.49012e-08, max_iter=10000):
    """
    Calculates the roots of the function through newton-raphson method.
    
    Parameters:
    -----------
    x0: initial guess
    N_x_tot: total number of moles in system
    t: time [yr]
    planet: planet
    star: star which planet orbits
    gas: species of gas
    tol: tolerance
    max_iter: number of the maximal iterations
        
    Returns:
    --------
    root : calculated root of the function
    """
    xold = x0 + 1
    xnew = x0
    numits = 0
    
    while (abs(xnew-xold) > tol):
        fxnew = func_mod.func(xnew, N_x_tot, planet, star, gas, melt_fraction)
        fprimenew = scipy.misc.derivative(func_mod.func, xnew, args = (N_x_tot, planet, star, gas, melt_fraction))
        numits = numits + 1
        
        xold = xnew
        xnew = xnew - fxnew/fprimenew
        
        root = xnew

    return root  # Returns the point at which function = 0


# Initialize classes
venus =  Planet(4.8675e24,6051800, 8.87, 0.7, 0.3, 1.07712e11, 2700, 620, 1e4)

# Gas properties of ideal gas at T = 300 K
co2 = Gas(188.92, 849, 0.0440095, 1.937e-9, 0.714)

sun = Star(3.828e26, 4.499e9)

# Conversion constants
bar_to_pascal = 1e5

N_x_tot = 1e22  # [moles]
melt_fraction_0 = 0.001  # [kg/kg]

# Plot f(pCO2) to look for roots 
x = np.logspace(-3, 8, 1000) # pressure [Pa]
f = func_mod.func(x, N_x_tot, venus, sun, co2, melt_fraction_0) # moles CO2 [mol]
# plt.loglog(x, f*(f>0),x, -f*(f<0))
# plt.title('Unique Species Partial Pressure for N_tot = {}'.format(N_x_tot))
# plt.xlabel('Pressure [Pa]', size='x-large')
# plt.ylabel('f($p_{H2O}$)', size='x-large')
# plt.grid()

x0 = newton(1e7, N_x_tot, venus, sun, co2, melt_fraction_0)
print('Root at:', x0)

# plt.scatter(x0, N_x_tot, label= 'pressure of $CO_2$ = {} Pa'.format(x0))
# plt.legend()
# plt.show()

# Now vary melt fraction and plot result
melt_fraction_a = np.logspace(-4,0,100)
pCO2_a = np.zeros([100,1])

for ii in range(0, 100):
    pCO2_a[ii] = newton(1e7, N_x_tot, venus, sun, co2, melt_fraction_a[ii])


# plt.semilogx(melt_fraction_a,pCO2_a)
# plt.xlabel('melt fraction [kg/kg]')
# plt.ylabel('p_{CO2} [Pa]')


sys.exit()







