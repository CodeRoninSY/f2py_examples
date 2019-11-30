#!/usr/bin/env python
"""
test.py
<2019-11-30> CodeRoninSY
Pass func as argument to fortran compiled module
$> f2py -c -m test test.f90
$> python test.py
"""
from __future__ import print_function
import math
import numpy as np
import test


if __name__ == "__main__":
    def func(x: float, t: float) -> float:
        return x * math.sin(t)

    x = [1.0, 2.0, 5.0]
    t = np.linspace(0, 2 * math.pi, 12)
    # z = 0.0
    z = test.test(x, t, func)

    print("z: {}".format(z))