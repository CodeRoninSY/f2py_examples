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

    print(f"size(t): {np.size(t)}, size(y): {np.size(y)}")

    z_eu = euler_v3.ode.euler(y0, t, f)
    print(f"z_eu = {z_eu}")

    z_he = euler_v3.ode.heun(y0, t, f)
    print(f"z_he = {z_he}")

    z_rk2a = euler_v3.ode.rk2a(y0, t, f)
    print(f"z_rk2a = {z_rk2a}")

    z_rk2b = euler_v3.ode.rk2b(y0, t, f)
    print(f"z_rk2b = {z_rk2b}")

    z_rk3 = euler_v3.ode.rk3(y0, t, f)
    print(f"z_rk3 = {z_rk3}")

    z_rk4_38 = euler_v3.ode.rk4_38(y0, t, f)
    print(f"z_rk4_38 = {z_rk4_38}")

    z_rk4 = euler_v3.ode.rk4(y0, t, f)
    print(f"z_rk4 = {z_rk4}")

    z_rk45, e = euler_v3.ode.rk45(y0, t, f)
    print(f"z_rk45 = {z_rk45}, e = {e}")
