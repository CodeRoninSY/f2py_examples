#!/usr/bin/env python
# euler_v1.py
#
# <2019-11-09> CodeRoninSY
# Demonstration example for Python calling Fortran90
# Compile euler_v1.f90 and import
from __future__ import print_function
import numpy as np
import math
from pprint import pprint
import matplotlib.pyplot as plt
import euler_v3

# print(euler_v3.ode.euler.__doc__)
# print(euler_v3.ode.heun.__doc__)


def f(y, t):
    return y * np.sin(t)


if __name__ == "__main__":
    a, b = (0.0, 10.0)
    y0 = -1.0
    n = 31
    tol = 1E-6
    hmax = 1.0
    hmin = 1E-6

    # initialize arrays
    t = np.linspace(a, b, n)
    y = np.array([y0]*n)
    y[0] = y0

    # analytical solution
    y_ex = -np.exp(1.0 - np.cos(t))

    print(f"size(t): {np.size(t)}, size(y): {np.size(y)}")

    y_eu = euler_v3.ode.euler(y0, t, f)
    print(f"y_eu = {y_eu}")

    y_he = euler_v3.ode.heun(y0, t, f)
    print(f"y_he = {y_he}")

    y_rk2a = euler_v3.ode.rk2a(y0, t, f)
    print(f"y_rk2a = {y_rk2a}")

    y_rk2b = euler_v3.ode.rk2b(y0, t, f)
    print(f"y_rk2b = {y_rk2b}")

    y_rk3 = euler_v3.ode.rk3(y0, t, f)
    print(f"y_rk3 = {y_rk3}")

    y_rk4_38 = euler_v3.ode.rk4_38(y0, t, f)
    print(f"y_rk4_38 = {y_rk4_38}")

    y_rk4 = euler_v3.ode.rk4(y0, t, f)
    print(f"y_rk4 = {y_rk4}")

    y_pc4 = euler_v3.ode.pc4(y0, t, f)
    print(f"y_pc4 = {y_pc4}")

    y_rk45, e = euler_v3.ode.rk45(y0, t, f)
    print(f"y_rk45 = {y_rk45}, e = {e}")

    #  figure( 1 )
    plt.figure(2, figsize=(14, 10))
    plt.subplot(2, 2, 1)
    plt.plot(t, y_ex, 'k-', t, y_eu, 'b--',
             t, y_he, 'g-*',
             t, y_rk2a, 'r-.', t, y_rk2b, 'y-.',
             t, y_rk3, 'r:')
    plt.xlabel('$t$')
    plt.ylabel('$y$')
    plt.title('Solutions of $dy/dt = y \sin t$, $y(0)=-1$')
    plt.legend(('Analytic', 'Euler', 'Heun', 'RK2A',
                'RK2B', 'RK3'),
               loc='lower left',
               framealpha=0.2, fontsize='small')

    # figure( 2 )
    plt.subplot(2, 2, 2)
    plt.plot(t, y_he - y_ex, 'g-*',
             t, y_rk2a - y_ex, 'r-.', t, y_rk2b - y_ex, 'y-.',
             t, y_rk3 - y_ex, 'r:')
    plt.xlabel('$t$')
    plt.ylabel('$y - y^*$')
    plt.title('Errors in solutions of $dy/dt = y \sin t$, $y(0)=-1$')
    plt.legend(('Heun', 'RK2A',
                'RK2B', 'RK3'),
               loc='upper left',
               framealpha=0.2, fontsize='small')

    # figure( 3 )
    plt.subplot(2, 2, 3)
    plt.plot(t, y_rk4, 'b-.', t, y_rk3, 'r-.',
             t, y_rk4_38,'k:' , t, y_rk45, 'g-.')
    plt.xlabel('$t$')
    plt.ylabel('$y$')
    plt.title('Solutions of $dy/dt = y \sin t$, $y(0)=-1$')
    plt.legend(('$O(h^4)$ Runge-Kutta','RK3', 'RK4-3/8',
                'RK45'),
               loc='lower left',
               framealpha=0.2, fontsize='small')

    # figure( 4 )
    plt.subplot(2, 2, 4)
    plt.plot(t, y_rk4 - y_ex, 'b-.', t, y_rk3 - y_ex, 'r-.',
             t, y_rk4_38 - y_ex, 'k:',t, y_rk45 - y_ex, 'g-.',
            )
    plt.xlabel('$t$')
    plt.ylabel('$y - y^*$')
    plt.title('Errors in solutions of $dy/dt = y \sin t$, $y(0)=-1$')
    plt.legend(('$O(h^4)$ Runge-Kutta','RK3', 'RK4-3/8',
                'RK45'),
               loc='upper left',
               framealpha=0.2, fontsize='small')

    # show()
    plt.savefig('euler_v3_Fig1.png', dpi=300)
