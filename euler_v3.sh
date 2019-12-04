#!/bin/bash

python -m numpy.f2py -c -m euler_v3 euler_v3.f90
python euler_v3.py
