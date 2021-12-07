'''
Created Thu Sep 30 2021
Collin Cherubim
'''

import numpy as np
import matplotlib as plt
plt.rcParams["figure.figsize"] = (8,6)

# Constants

T = 250  # temp in [K]
inv_cm2m = 100 # convert inverse cm to inverse m
avogadro = 6.022e23 # Avogadro's number [particles/mole]
s2day = 1/(3600*24) # convert s to day
s2yr = 1/(3600*24*365) # convert s to year
R = 8.314 # gas constant [J/mol/K]
M_N2 = 0.028014 # molar mass N2 [kg/mol]
M_H2 = 0.002016 # molar mass H2 [kg/mol]
# tau = 1019720994904 # [s]

# Set variable values

b = 2.7e17 * T**0.75 * inv_cm2m # binary diffusion coeff [molecules/m/s]
g = 9.8 # grav field strength [m/s/s]
m_bar = 0.02896/avogadro # average molecular mass [kg/molecule]
p = 1e5 # atmospheric surface pressure [Pa]
H_N2 = R*T/(M_N2*g) # N2 scale height [m]
H_H2 = R*T/(M_H2*g) # H2 scale height [m]

# Analytic solution

# Solve for e-folding time, tau
def tau(p, b, g, m_bar, H1, H2):
    return p/((b*g*m_bar)*(H2**-1 - H1**-1))

def H2_escape(t,N0):
    return N0*np.exp(-t/tau) # N0 is initial value

# Numerical solver - Forward Euler method

N0 = 1.17e49 # initial number of H2 molecules
tau = tau(p, b, g, m_bar, H_H2, H_N2) # e-folding time
delta_t = 0.001*tau # timestep
y0 = N0 # initialize

def f(y):
    return -y/tau

n_tot = 1000
y = y0

t_a = delta_t*np.linspace(0,n_tot,n_tot)
y_a = np.zeros(n_tot)

for n in range(n_tot):
    y = y + f(y) * delta_t
    y_a[n] = y

# Plot data

plt.plot(t_a*s2yr/1e6, np.log10(y_a), 'o', color = 'pink', label = 'numeric solution')
plt.plot(t_a*s2yr/1e6, np.log10(H2_escape(t_a,N0)), '--', label = 'analytic solution')
plt.title('Log(H2) vs Time', fontsize = 14)
plt.xlabel('Time [Ma]', fontsize = 14)
plt.ylabel('log($N_{H2}$) [molecules]', fontsize = 14)
plt.legend(fontsize = 14)
