#!/bin/bash

python -m numpy.f2py -c -m euler euler_v1.f90
python euler_v1.py
