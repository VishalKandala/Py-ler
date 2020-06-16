#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 10 17:55:48 2020
@author: vishal
"""
import numpy as np
# Solver Parameters #
d2r=np.pi/90 # 1 Degree to Radians
Imax=51 # No of grid points along x
Jmax=51 # No of grid points along Y
Icmax=Imax-1 #No of cells along Y
Jcmax=Jmax-1 # No of cells along X
Itmax=500 # Total number of Iterations.
nu2=0.5 # Second order dissipation constant
nu4=0.01 # Fourth order dissipation constant
cfl=0.5 # CFL Number
dt=0.1 # Global Timestep#
eps=1e-5 # Tolerance
# RK-4 Constants#
a=[0.25,0.33,0.5,1]
