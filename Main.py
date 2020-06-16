#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 10 13:48:06 2020
@author: vishal
"""
import numpy as np
import func as fc
import gasprop as gp
import flowprop as fp
import solverparam as slp
import flux as flx
import meshgen
import diffusion as df
###############################################################################
# Pre-Processing #
#mesh=meshgen.gen(slp.Imax,slp.Jmax,1000,'ps')
#X=mesh[0]
#Y=mesh[1]
# Read Mesh file #
Xmeshfile='X_mesh_51.txt'
Ymeshfile='Y_mesh_51.txt'
X=np.loadtxt(Xmeshfile)
Y=np.loadtxt(Ymeshfile)
# Cell Area/Volume
vol=np.zeros((slp.Icmax,slp.Jcmax))
for j in range(slp.Icmax):
for i in range(slp.Jcmax):
vol[i,j]=0.5*np.abs(((X[i+1,j+1]-X[i,j])*(Y[i,j+1]-Y[i+1,j]))-((Y[i+1,j+1]-Y[i,j])*(X[i,j+
# Connection between cells and nodes#
xnode=np.zeros((slp.Icmax,slp.Jcmax,4))
ynode=np.zeros((slp.Icmax,slp.Jcmax,4))
for i in range(slp.Icmax-1):
for j in range(slp.Jcmax):
xnode[i,j,0]=X[i,j]
ynode[i,j,0]=Y[i,j]
xnode[i,j,1]=X[i+1,j]
ynode[i,j,1]=Y[i+1,j]
xnode[i,j,2]=X[i+1,j+1]
ynode[i,j,2]=Y[i+1,j+1]
xnode[i,j,3]=X[i,j+1]
ynode[i,j,3]=Y[i,j+1]
# Faces dx and dy
xface=np.zeros((slp.Icmax,slp.Jcmax,4))
yface=np.zeros((slp.Icmax,slp.Jcmax,4))
for i in range(slp.Imax-1):
for j in range(slp.Jmax-1):
xface[i,j,0]=xnode[i,j,1]-xnode[i,j,0]
yface[i,j,0]=ynode[i,j,1]-ynode[i,j,0]
xface[i,j,1]=xnode[i,j,1]-xnode[i,j,2]
yface[i,j,1]=ynode[i,j,2]-ynode[i,j,1]
xface[i,j,2]=xnode[i,j,3]-xnode[i,j,2]
yface[i,j,2]=ynode[i,j,2]-ynode[i,j,3]
xface[i,j,3]=xnode[i,j,3]-xnode[i,j,0]
yface[i,j,3]=ynode[i,j,0]-ynode[i,j,3]
# Outward normal-Hence dx<0
# Outward normal dx<0
# Outward normal Hence dy<0
# Outward normal Hence dy<0
# Combined face vector
face=np.zeros((slp.Icmax,slp.Jcmax,4,2))
for k in range(4):
face[:,:,k,0]=xface[:,:,k]
face[:,:,k,1]=yface[:,:,k]
# Normals to faces
1fnormal=np.zeros((slp.Icmax,slp.Jcmax,4,2))
for k in range(4):
fnormal[:,:,k,0]=face[:,:,k,1]
fnormal[:,:,k,1]=(-1)*face[:,:,k,0]
facemag=np.zeros((slp.Icmax,slp.Jcmax,4))
for k in range(4):
facemag[:,:,k]=((face[:,:,k,0])**2+(face[:,:,k,1])**2)
# Local Eigen values #
# Local Maximum time step #
###############################################################################
# Variables # # Fill in Initial onditions.
# Primitive Variables#
rho=(fp.rinf)*np.ones((slp.Icmax,slp.Jcmax))
u=(fp.uinf)*np.ones((slp.Icmax,slp.Jcmax))
v=(fp.vinf)*np.ones((slp.Icmax,slp.Jcmax))
U=(fp.qinf)*np.ones((slp.Icmax,slp.Jcmax))
p=(fp.pinf)*np.ones((slp.Icmax,slp.Jcmax))
# For Eigenvalue #
c=(fp.cinf)*np.ones((slp.Icmax,slp.Jcmax))
u_abs=np.zeros((slp.Icmax,slp.Jcmax,4))
eigav=np.zeros((slp.Icmax,slp.Jcmax,4))
# State Vector #
q=np.zeros((slp.Icmax,slp.Jcmax,4))
q[:,:,0]=rho
q[:,:,1]=rho*u
q[:,:,2]=rho*v
q[:,:,3]=((p*gp.gamma_1)+((0.5*rho)*((u*u)+(v*v))))
# main pseudo-time loop #
###############################################################################
###############################################################################
res=[]
n=0
while n<(slp.Itmax):
dtmax=(slp.dt*np.ones((slp.Icmax,slp.Jcmax)))
# Impose Wall BCs
# Top and bottom boundaries #
#
v[:,slp.Jcmax-1]=0
#
q[:,slp.Jcmax-1,1]=q[:,slp.Jcmax-1,2]*(fnormal[:,slp.Jcmax-1,2,0]/fnormal[:,slp.Jcmax-1,2,1])
#
q[:,0,1]=q[:,0,2]*(fnormal[:,0,0,0]/fnormal[:,0,0,1])
#
v[:,0]=0
#
q[:,0,2]=0
p[:,0]=p[:,1]
p[:,slp.Jcmax-1]=p[:,slp.Jcmax-2]
# Flux BC's Implemented in fluxes module
###############################################################################
# Impose Inlet & Outlet Bcs
# Inlet
p[0,:]=p[1,:]
q[0,:,1]=rinf*uinf
q[0,:,2]=rinf*vinf
# Outlet
p[slp.Jcmax-1,:]=pinf
q[slp.Jcmax-1,:,1]=q[slp.Jcmax-2,:,1]
q[slp.Jcmax-1,:,2]=q[slp.Jcmax-2,:,2]
q[slp.Jcmax-1,:,3]=q[slp.Jcmax-2,:,3]
###############################################################################
D=df.diff(q,face,fnormal)# Calculate Damping
###############################################################################
# Runge-Kutta loop #
###############################################################################
qrk=np.zeros((5,slp.Icmax,slp.Jcmax,4))
qrk[0]=q
rrk=np.zeros((5,slp.Icmax,slp.Jcmax,4))
rrk[0]=flx.R(rho,u,v,U,p,q,fnormal)
for k in range(4):
u_abs[:,:,k]=np.abs((u[:,:]*fnormal[:,:,k,0])+(v[:,:]*fnormal[:,:,k,1]))
c[:,:]=((gp.gamma)*(p[:,:]/rho[:,:]))
2cface=np.zeros((slp.Icmax,slp.Jcmax,4))
for i in range(1,slp.Icmax-1):
for j in range(1,slp.Jcmax-1):
cface[i,j,0]=0.5*(c[i,j]+c[i,j-1])
cface[i,j,1]=0.5*(c[i,j]+c[i+1,j])
cface[i,j,2]=0.5*(c[i,j+1]+c[i,j])
cface[i,j,3]=0.5*(c[i,j]+c[i-1,j])
for k in range(4):
eigav[:,:,k]=u_abs[:,:,k]+cface[:,:,k]
for k in range(1,5):
for i in range(1,slp.Icmax-1):
for j in range(1,slp.Jcmax-1):
den=np.ones((slp.Icmax-1,slp.Jcmax))
for l in range(4):
den[i,j]+= np.abs(eigav[i,j,k-1])*facemag[i,j,k-1]
dtmax[i,j]=2*vol[i,j]/(den[i,j])
qrk[k,i,j]=qrk[k-1,i,j]-(slp.a[k-1]*((slp.cfl*dtmax[i,j])/vol[i,j])*(rrk[k-1,i,j])
rho[:,:]=qrk[k,:,:,0]
u[:,:]=qrk[k,:,:,1]/rho[:,:]
v[:,:]=qrk[k,:,:,2]/rho[:,:]
U[:,:]=((u[:,:])**2+(v[:,:])**2)**0.5
p[:,:]=gp.gamma1*(qrk[k,:,:,3]-((1/(2*qrk[k,:,:,0]))*((qrk[k,:,:,1]**2)+(qrk[k,:,:,2]**2))
rrk[k]=flx.R(rho,u,v,U,p,qrk[k],fnormal)
R=rrk[k]
###############################################################################
# Residual #
resvec=np.zeros((slp.Icmax,slp.Jcmax,4))
resvec[:,:]=np.abs(D[:,:]-R[:,:])
resmax=np.amax(resvec)
iterres=np.array([n,resmax])
res.append(iterres)
#Check convergence #
err=np.zeros((slp.Icmax,slp.Jcmax,4))
err=qrk[4]-qrk[0]
ermax=np.amax(err)
if ermax<slp.eps:
break
###############################################################################
q=qrk[4]
n+=1
3
