#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 11 02:48:14 2020
@author: vishal
"""
import numpy as np
import solverparam as slp
import flowprop as fp
import gasprop as gp
def switch2(p):
s2=np.zeros((2,slp.Icmax,slp.Jcmax))
for i in range(slp.Icmax-1):
for j in range(slp.Jcmax-1):
if j==slp.Jcmax:
pden=(1/(p[i+1,j]+(2*p[i,j])+p[i-1,j]))
pnum1=(p[i,j]-(2*p[i,j])+p[i,j-1])
pnum0=(p[i+1,j]-(2*p[i,j])+p[i-1,j])
elif j==0:
pden=(1/(p[i+1,j]+(2*p[i,j])+p[i-1,j]))
pnum0=(p[i+1,j]-(2*p[i,j])+p[i-1,j])
pnum1=(p[i,j+1]-(2*p[i,j])+p[i,j])
else:
pden=(1/(p[i+1,j]+(2*p[i,j])+p[i-1,j]))
pnum0=(p[i+1,j]-(2*p[i,j])+p[i-1,j])
pnum1=(p[i,j+1]-(2*p[i,j])+p[i,j-1])
s2[0,i,j]=slp.nu2*abs(pnum0)*pden
s2[1,i,j]=slp.nu2*abs(pnum1)*pden
s2face=np.zeros((2,4,slp.Icmax,slp.Jcmax))
for i in range(slp.Icmax-1):
for j in range(slp.Jcmax-1):
s2face[0,0,i,j]=0.5*(s2[0,i,j]+s2[0,i,j-1])
s2face[0,1,i,j]=0.5*(s2[0,i,j]+s2[0,i+1,j])
s2face[0,2,i,j]=0.5*(s2[0,i,j+1]+s2[0,i,j])
s2face[0,3,i,j]=0.5*(s2[0,i,j]+s2[0,i-1,j])
##########################################
s2face[1,0,i,j]=0.5*(s2[1,i,j]+s2[1,i,j-1])
s2face[1,1,i,j]=0.5*(s2[1,i,j]+s2[1,i+1,j])
s2face[1,2,i,j]=0.5*(s2[1,i,j+1]+s2[1,i,j])
s2face[1,3,i,j]=0.5*(s2[1,i,j]+s2[1,i-1,j])
return s2face
def diff(q,face,fnormal):
D=np.zeros((slp.Icmax,slp.Jcmax,4))
Dface=np.zeros((4,slp.Icmax,slp.Jcmax,4))
facemag=np.zeros((slp.Icmax,slp.Jcmax,4))
u=(fp.uinf)*np.ones((slp.Icmax,slp.Jcmax))
v=(fp.vinf)*np.ones((slp.Icmax,slp.Jcmax))
p=(fp.pinf)*np.ones((slp.Icmax,slp.Jcmax))
c=(fp.cinf)*np.ones((slp.Icmax,slp.Jcmax))
u_abs=np.zeros((slp.Icmax,slp.Jcmax,4))
eigav=np.zeros((slp.Icmax,slp.Jcmax,4))
rho=(fp.rinf)*np.ones((slp.Icmax,slp.Jcmax))
rho[:,:]=q[:,:,0]
u[:,:]=q[:,:,1]/rho[:,:]
v[:,:]=q[:,:,2]/rho[:,:]
p[:,:]=gp.gamma1*(q[:,:,3]-((1/(2*q[:,:,0]))*((q[:,:,1]**2)+(q[:,:,2]**2))))
for k in range(4):
facemag[:,:,k]=((face[:,:,k,0])**2+(face[:,:,k,1])**2)
u_abs[:,:,k]=np.abs((u[:,:]*fnormal[:,:,k,0])+(v[:,:]*fnormal[:,:,k,1]))
c[:,:]=((gp.gamma)*(p[:,:]/rho[:,:]))
cface=np.zeros((slp.Icmax,slp.Jcmax,4))
for i in range(slp.Icmax-1):
for j in range(slp.Jcmax-1):
cface[i,j,0]=0.5*(c[i,j]+c[i,j-1])
1cface[i,j,1]=0.5*(c[i,j]+c[i+1,j])
cface[i,j,2]=0.5*(c[i,j+1]+c[i,j])
cface[i,j,3]=0.5*(c[i,j]+c[i-1,j])
for k in range(4):
eigav[:,:,k]=u_abs[:,:,k]+cface[:,:,k]
s2=switch2(p)
s4=np.zeros((2,4,slp.Icmax,slp.Jcmax))
s4[0,:,:]=np.maximum(0,(slp.nu4-s2[0,:,:]))
dq=np.zeros((4,slp.Icmax,slp.Jcmax,4))
d3q=np.zeros((4,slp.Icmax,slp.Jcmax,4))
for i in range(1,slp.Icmax-2):
for j in range(1,slp.Jcmax-2):
dq[0,i,j]=q[i,j]-q[i,j-1]
dq[2,i,j]=q[i,j+1]-q[i,j]
dq[1,i,j]=q[i+1,j]-q[i,j]
dq[3,i,j]=q[i,j]-q[i-1,j]
d3q[2,i,j]=q[i,j+2]-2*q[i,j+1]+q[i,j]-(q[i,j+1]-2*q[i,j]+q[i,j-1])
d3q[1,i,j]=q[i+2,j]-2*q[i+1,j]+q[i,j]-(q[i+1,j]-2*q[i,j]+q[i-1,j])
d3q[0,i,j]=q[i,j+1]-2*q[i,j]+q[i,j-1]-(q[i,j]-2*q[i,j-1]+q[i,j-2])
d3q[3,i,j]=q[i+1,j]-2*q[i,j]+q[i-1,j]-(q[i,j]-2*q[i-1,j]+q[i-2,j])
for i in range(slp.Icmax):
for j in range(slp.Jcmax):
for k in range(4):
Dface[k,i,j]+=(s2[0,k,i,j]*eigav[i,j,k]*facemag[i,j,k]*dq[k,i,j])-(s4[0,k,i,j]*eig
D[:,:]=Dface[0,:,:]+Dface[1,:,:]+Dface[2,:,:]+Dface[3,:,:]
return D
2
