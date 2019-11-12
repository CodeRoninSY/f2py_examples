#!/usr/bin/env python
# newt_raph.py
#
# <2019-11-09> CodeRoninSY
# Demonstration example for Python calling Fortran90
# Compile newt_raph.f90 and import
from __future__ import print_function
import numpy as np
import math
import newt_raph as nr

# print(nr.newton.solve.__doc__)

def f_sqrt(x):
    return x ** 2 - 4.0

def fprime_sqrt(x):
    return 2.0 * x

def f1(x):
    return math.cos(x) + 2 * math.sin(x) + x * x

def f1pr(x):
    return -1 * math.sin(x) + 2 * math.cos(x) + 2 * x

def f2(x):
    return math.exp(-x) * math.cos(x)

def f2pr(x):
    return - 1 * math.exp(-x) * (math.sin(x) + math.cos(x))


def main():
    x = 0       # (out) x
    iters = 0   # (out) iters
    debug = True

    print(f"Test routine for computing zero of f")
    x0vals = [-1.5, -1.0, -0.6598, -0.5, -0.1, 0.1, 1.0, 1.3]

    for itest in range(len(x0vals)):
        x0 = x0vals[itest]
        x, iters = nr.newton.solve(f2, f2pr, x0, x, iters, debug)
        print(f"solver returns x= {x:22.12e} after {iters:4d}")
        fx = f2(x)
        print(f"the value of f(x) is {fx}")
        # if (abs(x - 2.0) > 1.0e-6):
        #     print(f"*** unexpected result: x = {x}")


if __name__ == "__main__":
    main()
