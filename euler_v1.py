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
import euler

# print(euler.ode.euler.__doc__)
# print(euler.ode.heun.__doc__)

# function kwargs for solvers
paramsNWCOL = {
    "y0": 100.0,
    "a": 0.0,
    "b": 200.0,
    "h": 10.0
}

params = {
    "y0": -1.0,
    "a": 0.0,
    "b": 10.0,
    "h": .1
}

# function kwargs for RKF
paramRKF = {
    "y0": -1.0,
    "a": 0.0,
    "b": 10.0,
    "tol": 1E-6,
    "hmax": 1.0,
    "hmin": 1E-2
}


def f(x, t):
    return x * np.sin(t)


if __name__ == "__main__":
    a, b = (0.0, 10.0)
    y0 = -1.0
    n = 101
    tol = 1E-6
    hmax = 1.0
    hmin = 1E-6
    
    # initialize arrays
    t = np.linspace(a, b, n)
    y = np.array([y0]*n)
    y[0] = y0

    print(f"size(t): {np.size(t)}, size(y): {np.size(y)}")

    z_eu = euler.ode.euler(y0, t, f)
    print(f"z_eu = {z_eu}")

    z_he = euler.ode.heun(y0, t, f)
    print(f"z_he = {z_he}")

    z_rk4 = euler.ode.rk4(y0, t, f)
    print(f"z_rk4 = {z_rk4}")

    z_rk45, e = euler.ode.rk45(y0, t, f)
    print(f"z_rk45 = {z_rk45}, e = {e}")

    z_rkf, t_ = euler.ode.rkf(y0, a, b, tol, hmin, hmax, f)
    print(f"z_rkf = {z_rkf}, t_ = {t_}")
