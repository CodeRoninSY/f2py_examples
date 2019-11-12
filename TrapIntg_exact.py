#!/usr/bin/env python
from __future__ import print_function
import numpy as np
import math as m
import matplotlib.pyplot as plt
import time

# THE EXACT SOLUTION
# OBS: Takes x/L as argument, instead of only x

Phi_Res = 1000 # Resolution in phi: no. of times to perform the integral
Int_Res = 1000 # Resolution in the integral: Trapezoidal method

L = 1.0

x_list   = np.linspace(0, L, Phi_Res)
phi_x_ex = np.zeros(Phi_Res)


def ExactSolution(xL):
    Term1 = m.sqrt((xL - L)**2 + 1) - xL*m.asinh(xL - L) # z = L
    Term2 = m.sqrt((xL - 0)**2 + 1) - xL*m.asinh(xL - 0) # z = 0
    return Term1 - Term2

for i in range(Phi_Res):
    phi_x_ex[i] = ExactSolution(x_list[i] / L)

# CONSTANTS
PI      = np.pi     # radians [dimless]
L       = 1.0       # meters
lambda0 = 0.3       # Coloumb/meters
eps0 = 8.854e-11  # Farads/meters

phi_x    = np.zeros(Phi_Res) # initialize to zero

z_list   = np.linspace(0, L, Int_Res)

dz       = z_list[1] - z_list[0] # delta_z

def chargeDens(z_):
    # lambda= lambda0 * x' = (lambda0*L)*z
    # and lambda0*L is taken outside the integral
    # such that we're left with:
    return z_

def integrand(z_, x_):
    return (chargeDens(z_)) / (np.sqrt((x_/L - z_)**2 + 1))


start = time.time()
i = 0
for x in x_list: # Perform the integral for all x-values
    for z in z_list:
        phi_x[i] = phi_x[i] + (integrand(z + dz, x) + integrand(z, x))
    phi_x[i] = phi_x[i] * 0.5 * dz
    i = i + 1
end = time.time()
purePythonTime = end-start
print("Elapsed time (s): ", round(purePythonTime,3))


fontSettings = {'size' : 16 }
plt.rc('font', **fontSettings)

plt.figure(figsize=(10, 6))
plt.title('Exact vs (Pure Python) Numerical Solution')
plt.ylabel('C Psi(z)')
plt.xlabel('z=x/L')
plt.plot(x_list,phi_x_ex,color='g', linestyle='-', linewidth=3.0, label='Exact')
plt.plot(x_list,phi_x,   color='r', linestyle='--',linewidth=3.0, label='Numerical')
plt.legend(loc=4)
plt.show()


import TrapIntg_elRod as TIeR
# The Python module containing compiled Fortran code.
# The Fortran code is included at the very end of this Example.


# Check the documentation of out Fortran Module:
print("--------------------------------")
print(TIeR.__doc__)
print("--------------------------------")


start = time.time()
phi_x = TIeR.trapezoidalintegrator(np.array(Phi_Res), np.array(Int_Res), L)
end = time.time()

F2PYtime = end-start
print("Elapsed time (s): ", round(F2PYtime, 4))

plt.figure(figsize=(10, 6))
plt.title('Exact vs (F2PY) Numerical Solution')
plt.ylabel('C Psi(z)')
plt.xlabel('z=x/L')
plt.plot(x_list, phi_x_ex, color='g', linestyle='-',
         linewidth=3.0, label='Exact')
plt.plot(x_list, phi_x,   color='r', linestyle='--',
         linewidth=3.0, label='Numerical')
plt.legend(loc=4)
plt.show()

print("PurePythonTime / F2PYtime ~= ", round(purePythonTime/F2PYtime))
print("\nThe elapsed F2PY-time constitutes only", round((F2PYtime/purePythonTime)*100, 2),
      "% of the elapsed pure Python time.")
