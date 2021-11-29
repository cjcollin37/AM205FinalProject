#!/usr/bin/env python3
# -*- coding: utf-8 -*-
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






# Modeling Magma Ocean and Atmospheric Interaction on Venus
# Import relevant libraries
import numpy as np
import matplotlib.pyplot as plt
import scipy.misc


# Import modules
import func
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





# Initialize classes
venus =  Planet(4.8675e24,6051800, 8.87, 0.7, 0.3, 1.07712e11, 2700, 620, 1e4)

# Gas properties of ideal gas at T = 300 K
co2 = Gas(188.92, 849, .0440095, 1.937e-9, 0.714)
co2_exp = Gas(188.92, 849, .0440095, 1.25e-3, 0.714) # CO2 data with experimental alpha value (expect partitioning to be roughly equal at this alpha)
h2o = Gas(461.5, 1872.3, .018015, 1.033e0, 1.747)
h2 = Gas(412.4, 14307, .002016, 2.572e-6, 1.000)
ch4 = Gas(518.2, 2253.7, .016043, 9.937e-8, 1.000)
n2 = Gas(296.8, 1039, .028013, 7.000e-5, 1.800)
co = Gas(296.8, 1040, .028011, 1.6e-7, 1.000)

sun = Star(3.828e26, 4.499e9)

# Conversion constants
bar_to_pascal = 1e5

# Initialize other variables
# Time array
t = np.linspace(0,4.4999e9) # [yr]
# Pressure array
p_s = np.linspace(0.01, 1e7) # [Pa]


N_x_tot = 1e22 #[moles]


# Newton's method function to solve for pressure
def newton(x0, N_x_tot, planet, star, gas, tol=1.49012e-08, max_iter=10000):
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
        fxnew = func.func(xnew, N_x_tot, planet, star, gas)
        fprimenew = scipy.misc.derivative(func.func, xnew, args = (N_x_tot, planet, star, gas))
        numits = numits + 1
        
        xold = xnew
        xnew = xnew - fxnew/fprimenew
        
        root = xnew

    return root  # Returns the point at which function = 0






# # * Plotting f(pCO2) to look for roots (=0)
x = np.linspace(0.001, 1e12) # pressure [Pa]
f = func.func(x, N_x_tot, venus, sun, h2)
plt.plot(x, f)
plt.title('Unique Species Partial Pressure for N_tot = {}'.format(N_x_tot))
plt.xlabel('Pressure [Pa]', size='x-large')
plt.ylabel('f($p_{H2O}$)', size='x-large')
plt.grid()


# root = newton(1e7, N_x_tot, venus, sun, co2)
root = newton(1e7, N_x_tot, venus, sun, co2_exp)
print('Root at:', root)

plt.scatter(root, 0, label= 'pressure of $CO_2$ = {} Pa'.format(root))
plt.legend()
plt.show()

N_atm_root, N_mo_root = moles.partitioning(root, venus, sun, co2_exp)
print('Root returns:', N_atm_root, N_mo_root)









# * Calculate partitioning from moles function
N_atm, N_mo = moles.partitioning(p_s, venus, sun, h2o)
Total = N_atm + N_mo
print('Total = ', Total)

# * Plotting Partitioning by surface pressure
plt.plot(p_s, N_atm, label='Moles in atmosphere')
plt.plot(p_s, N_mo, label='Moles in magma ocean')
plt.title('Partitioning vs. Surface Pressure')
plt.xlabel('$P_{surf}$ [Pa]')
plt.ylabel('$CO_2$ Abundance [moles]')
plt.legend()
plt.grid()
plt.show()






# # * Plotting  Magma Ocean Partitioning by pressure for different species
N_atm_co2, N_mo_co2 = moles.partitioning(p_s, venus, sun, co2)
N_atm_h2o, N_mo_h2o = moles.partitioning(p_s, venus, sun, h2o)
N_atm_h2, N_mo_h2 = moles.partitioning(p_s, venus, sun, h2)
N_atm_ch4, N_mo_ch4 = moles.partitioning(p_s, venus, sun, ch4)
N_atm_n2, N_mo_n2 = moles.partitioning(p_s, venus, sun, n2)
N_atm_co, N_mo_co = moles.partitioning(p_s, venus, sun, co)


total_co2 = N_atm_co2 + N_mo_co2
total_h2o = N_atm_h2o + N_mo_h2o
total_co = N_atm_co + N_mo_co
total_h2 = N_atm_h2 + N_mo_h2
total_ch4 = N_atm_ch4 + N_mo_ch4
total_n2 = N_atm_n2 + N_mo_n2

plt.plot(p_s, N_mo_co2, label='CO2', color='firebrick')
# plt.plot(p_s, N_mo_h2o, label='H2O', color='steelblue')
plt.plot(p_s, N_mo_h2, label='H2', color='green')
plt.plot(p_s, N_mo_ch4, label='CH4', color='mediumvioletred')
plt.plot(p_s, N_mo_n2, label='N2')
plt.plot(p_s, N_mo_co, label='CO', color='peru')
plt.title('MO Partitioning vs. Pressure')
plt.xlabel('Pressure [Pa]')
plt.ylabel('Abundance of species in magma ocean[moles]')
plt.legend()
plt.grid()
plt.show()




## Lol is any of this correct
# * Plotting Total Partitioning over time for different species
# plt.loglog(t, N_mo_co2, label='MO CO2', color='firebrick')
# plt.loglog(t, N_atm_co2, '--', label='Atm CO2', color='firebrick')
# plt.loglog(t, N_mo_h2o, '--', label='MO H2O', color='steelblue')
# plt.loglog(t, N_atm_h2o, '--', label='Atm H2O', color='steelblue')
# plt.loglog(t, N_mo_h2, label='MO H2', color='green')
# plt.loglog(t, N_atm_h2,'--', label='Atm H2', color='green')
# plt.loglog(t, N_mo_ch4, label='MO CH4', color='mediumvioletred')
# plt.loglog(t, N_mo_ch4, label='Atm CH4', color='mediumvioletred')
# plt.loglog(t, N_mo_n2, label='MO N2')
# plt.loglog(t, N_atm_n2, '--', label='Atm N2')
# plt.loglog(t, N_mo_co, label='MO CO', color='peru')
# plt.loglog(t, N_atm_co, '--', label='Atm CO', color='peru')
# plt.title('Partitioning over time')
# plt.xlabel('log(Time) [Yr]')
# plt.ylabel('Abundance [moles]')
# plt.legend(bbox_to_anchor=(1.3, 0.5), loc='right')
# plt.grid()
# plt.show()





# * From Lichtenberg et. al 2021 eq. 9 & Table 1
p_s = np.linspace(0.01, 1e6) # [bar]
x_i_co2 = co2.alpha*((p_s*bar_to_pascal)**(1/co2.beta))
x_i_h2o = h2o.alpha*((p_s*bar_to_pascal)**(1/h2o.beta))
x_i_h2 = h2.alpha*((p_s*bar_to_pascal)**(1/h2.beta))
x_i_ch4 = ch4.alpha*((p_s*bar_to_pascal)**(1/ch4.beta))
x_i_n2 = n2.alpha*((p_s*bar_to_pascal)**(1/n2.beta))
x_i_co = co.alpha*((p_s*bar_to_pascal)**(1/co.beta))

# * Recreate Lichtenberg et. al 2021 Fig.2
plt.loglog(p_s, x_i_h2o, label='H2O', color='steelblue')
plt.loglog(p_s, x_i_h2, label='H2', color='green')
plt.loglog(p_s, x_i_co2, label='CO2', color='firebrick')
plt.loglog(p_s, x_i_ch4, label='CH4', color='mediumvioletred')
#plt.loglog(p_s, x_i_n2, label='N2')
plt.loglog(p_s, x_i_co,'--', label='CO', color='peru')

# * Plot formatting
plt.ylim(1e-6,1e6)
plt.xlim(0.01,1e6)
plt.xlabel(' Gas Phase Pressure, $p_{vapor}^i [bar]$')
plt.ylabel('Volatile Solubility, $X_{magma}^i$ [ppmw]')
plt.legend()
