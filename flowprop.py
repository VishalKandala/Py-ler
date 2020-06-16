#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 10 13:34:08 2020
@author: vishal
"""
import solverparam as slp
import gasprop as gp
import numpy as np
M=0.3 # Mach Number @ inf
rinf=1.25# Rho @ inf
pinf=101325# Pressure @ inf
tinf=298 # static Temp @ inf
#ptinf= # Stagnation Pressure @inf
cinf=(gp.gamma*gp.Rc*tinf)**0.5 # Speed of sound @ inf
qinf=M*cinf# q @ inf
alpha=0 # Incidence @ inf
alpha_rad=alpha*slp.d2r
uinf=qinf*np.cos(alpha_rad) # u @ inf
vinf=qinf*np.sin(alpha_rad)# v @ inf
pratio=3 # Exit Static Pressure/ Inlet Stagnation Pressure
#pexit= # Static Pressure @ exit
#htinf= # Stagnation Enthalpy @ inf
#epsinf=# rho_E @ inf
1
