#!/usr/bin/env python
# euler_v1.py
#
# <2019-11-09> CodeRoninSY
# Demonstration example for Python calling Fortran90
# Compile euler_v1.f90 and import
from __future__ import print_function
import numpy as np
import math
import euler_v1

def newtoncooling(temp, t):
    return -0.0768 * (temp - 20)


def main():
    euler_v1.euler(newtoncooling,100, 0, 200, 5.00)


if __name__ == "__main__":
    main()
