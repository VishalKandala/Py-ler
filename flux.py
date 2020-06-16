#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 10 23:44:42 2020
@author: vishal
"""
import numpy as np
import solverparam as slp
import gasprop as gp
def R(rho,u,v,U,p,q,fnormal):
rho[:,:]=q[:,:,0]
u[:,:]=q[:,:,1]/rho[:,:]
v[:,:]=q[:,:,2]/rho[:,:]
U[:,:]=((u[:,:])**2+(v[:,:])**2)**0.5
p[:,:]=gp.gamma1*(q[:,:,3]-((1/(2*q[:,:,0]))*((q[:,:,1]**2)+(q[:,:,2]**2))))
# Flux variables and vectors#
f=np.zeros((slp.Icmax,slp.Jcmax,4))
g=np.zeros((slp.Icmax,slp.Jcmax,4))
#####
f[:,:,0]=q[:,:,1]
f[:,:,1]=(rho[:,:]*(u[:,:])**2)+p[:,:]
f[:,:,2]=(rho[:,:]*v[:,:])
f[:,:,3]=((q[:,:,3]*u[:,:])+(0.5*(U[:,:])**2/rho[:,:]))
f[:,0,0]=0
f[:,0,1]=p[:,0]*fnormal[:,0,0,1]
f[:,0,2]=0
f[:,0,3]=0
f[:,slp.Jcmax-1,0]=0
f[:,slp.Jcmax-1,1]=p[:,0]*fnormal[:,0,0,1]
f[:,slp.Jcmax-1,2]=0
f[:,slp.Jcmax-1,3]=0
########## Y flux ################
g[:,:,0]=q[:,:,2]
g[:,:,2]=(((q[:,:,2]**2)/q[:,:,0])+p[:,:])
g[:,:,1]=(rho[:,:]*u[:,:]*v[:,:])
g[:,:,3]=((q[:,:,3]*v[:,:])+(0.5*(q[:,:,2])**3/rho[:,:]**2))
g[:,0,0]=0
g[:,0,1]=0
g[:,0,2]=p[:,0]*fnormal[:,0,0,0]
g[:,0,3]=0
# Face fluxes #
fface=np.zeros((slp.Icmax,slp.Jcmax,4,4))
gface=np.zeros((slp.Icmax,slp.Jcmax,4,4))
for i in range(1,slp.Icmax-1):
for j in range(1,slp.Jcmax-1):
fface[i,j,0]=0.5*(f[i,j]+f[i,j-1])
fface[i,j,1]=0.5*(f[i,j]+f[i+1,j])
fface[i,j,2]=0.5*(f[i,j+1]+f[i,j])
fface[i,j,3]=0.5*(f[i,j]+f[i-1,j])
##################################
gface[i,j,0]=0.5*(g[i,j]+g[i,j-1])
gface[i,j,1]=0.5*(g[i,j]+g[i+1,j])
gface[i,j,2]=0.5*(g[i,j+1]+g[i,j])
gface[i,j,3]=0.5*(g[i,j]+g[i-1,j])
############################################################################
R=np.zeros((slp.Icmax,slp.Jcmax,4))
for i in range(slp.Icmax):
for j in range(slp.Jcmax):
for k in range(4):
R[i,j]+=((fface[i,j,k]*fnormal[i,j,k,1])+(gface[i,j,k]*fnormal[i,j,k,0]))
return R
1
