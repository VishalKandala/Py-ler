#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 10 14:03:24 2020
@author: vishal
"""
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 19:48:21 2020
@author: vishal
"""
import numpy as np
import matplotlib.pyplot as plt
#import math as m
#Functions#
def yl(x):
    n=len(x)
    yl=np.zeros(n)
    for i in range(n):
    if (x[i]>=2):
    if (x[i]<=3):
    yl[i]=0.1*np.sin((x[i]-2)*np.pi)
    return yl
def gen(nx,ny,T,p):
#
#
#
#########################
nx=51 # No of x-grid points(Imax)
ny=51 # No of y-grid points(Jmax)
T=2000 # No. of Iterations(Tmax)
dn=1 # Computational Domain x-difference #
de=1 # Computational Domain y-difference #
#erx=1 # Error in successive iterations X-coordinate
#ery=1 # Error in successive iterations Y-coordinate
#erc=1e-5 # Minimum error for convergence
#########################################
#Parameter Initialisation#
a=np.zeros((nx,ny)) # Alpha Matrix
b=np.zeros((nx,ny)) #Beta Matrix
g=np.zeros((nx,ny)) #Gamma Matrix
d=np.zeros((nx,ny)) # Denominator Parameter Matrix (Ease of Coding)
X0=np.zeros((nx,ny)) # X Co-ordinate Matrix
Y0=np.zeros((nx,ny)) # Y Co-ordinate Matrix
# Physical Domain Boundaries #
xmin=0
xmax=5
ymin=0
ymax=1
# Boundary Conditions #
c1=[np.linspace(xmin,xmax,nx),np.repeat(ymax,nx)] #Top Wall
c2=[np.repeat(xmin,ny),np.linspace(ymin,ymax,ny)]#Left Wall
c3=[np.repeat(xmax,ny),np.linspace(ymin,ymax,ny)]#Right Wall
c4=[np.linspace(xmin,xmax,nx),yl(np.linspace(xmin,xmax,nx))]#Bottom Wall
#
#Grid Initialisation#
X0[:,ny-1]=c1[0] #Top Wall
X0[0,:]=c2[0] #Left Wall
X0[nx-1,:]=c3[0] #Right Wall
X0[:,0]=c4[0] #Bottom Wall
Y0[:,ny-1]=c1[1] #Top Wall
Y0[0,:]=c2[1] #Left Wall
Y0[nx-1,:]=c3[1] #Right Wall
Y0[:,0]=c4[1] #Bottom Wall
1#Temporary Variables#
x=X0
y=Y0
Xt=x
Yt=y
#Iteration loop Begins#
for t in range(T): # Iterations loop
for i in range(1,nx-1): #X-axis sweep
for j in range(1,ny-1):#Y-axis sweep
# Parameter Computation (Frozen,Lagged)#
a[i,j]= (1/(4*(dn**2)))*(((x[i,j+1]-x[i,j-1])**2)+((y[i,j+1]-y[i,j-1])**2)) #Alpha
b[i,j]=(1/(4*de*dn))*(((x[i+1,j]-x[i-1,j])*(x[i,j+1]-x[i,j-1]))+((y[i+1,j]-y[i-1,j
g[i,j]=(1/(4*(de**2)))*(((x[i+1,j]-x[i-1,j])**2)+((y[i+1,j]-y[i-1,j])**2)) #Gamma
d[i,j]=((a[i,j]/(de**2))+(g[i,j]/(dn**2))) # Denominator for the following equatio
## Main Iteration of Laplace Equation #
Xt[i,j]=( (0.5/d[i,j])*( ( (a[i,j]/de**2)*(x[i+1,j]+x[i-1,j]) )-( (b[i,j]/(2*dn*de
Yt[i,j]=( (0.5/d[i,j])*( ( (a[i,j]/de**2)*(y[i+1,j]+y[i-1,j]) )-( (b[i,j]/(2*dn*de
#
Xt[:,ny-1]=c1[0] #Top Wall
#
Xt[0,:]=c2[0] #Left Wall
#
Xt[nx-1,:]=c3[0] #Right Wall
#
Xt[:,0]=c4[0] #Bottom Wall
#
Yt[:,ny-1]=c1[1] #Top Wall
#
Yt[0,:]=c2[1] #Left Wall
#
Yt[nx-1,:]=c3[1] #Right Wall
#
Yt[:,0]=c4[1] #Bottom Wall
# Neumann Boundary Conditions
if t>=(nx*ny/10): # Giving the equation some time to develop initially
Yt[1,1:ny-2]=Yt[0,1:ny-2] # Left Boundary
Yt[nx-2,1:ny-2]=Yt[nx-1,1:ny-2] # Right Boundary
Xt[1:nx-2,ny-2]=Xt[1:nx-2,ny-1]# Top Boundary
for n in range(1,nx-1): # Bottom Boundary Loop
if (Xt[n,0])>=2 and (Xt[n,0])<=2.25:
Xt[n,1]=Xt[n,0]#+(de*0.2*np.sin(np.pi*(Xt[n,0]-2)))
#
Yt[n,1]=Yt[n,0]+((1/(2*np.pi)*dn*0.2*np.cos(np.pi*(Xt[n,0]-2))))
elif (Xt[n,0])>2.25 and (Xt[n,0])<=2.5:
Xt[n,1]=Xt[n,0]#-(de*0.2*np.cos(np.pi*(Xt[n,0]-2)))
#
Yt[n,1]=Yt[n,0]+((1/(2*np.pi)*dn*0.2*np.sin(np.#pi*(Xt[n,0]-2))))
elif (Xt[n,0])>2.5 and (Xt[n,0])<=2.75:
Xt[n,1]=Xt[n,0]#+(de*0.2*np.cos(np.pi*(Xt[n,0]-2)))
#
Yt[n,1]=Yt[n,0]+((1/(2*np.pi)*dn*0.2*np.sin(np.pi*(Xt[n,0]-2))))
elif (Xt[n,0])>2.75 and (Xt[n,0])<=3:
Xt[n,1]=Xt[n,0]#+(de*0.2*np.sin(np.pi*(Xt[n,0]-2)))
#
Yt[n,1]=Yt[n,0]-((1/(2*np.pi)*dn*0.2*np.cos(np.pi*(Xt[n,0]-2))))
else:
Xt[n,1]=Xt[n,0] # X=(0,2) and (3,5) Flat boundaries.
## Variable Swap
x=Xt
y=Yt
if p=='p':
#Plot
r=plt.figure(1)
for i in range(nx):
plt.plot(Xt[i,:],Yt[i,:],'b')
for j in range(ny):
plt.plot(Xt[:,j],Yt[:,j],'g')
r.show()
elif p=='ps':
r=plt.figure(1)
for i in range(nx):
plt.plot(Xt[i,:],Yt[i,:],'b')
for j in range(ny):
plt.plot(Xt[:,j],Yt[:,j],'g')
r.show()
np.savetxt('X_mesh_'+str(nx)+'.txt',Xt)
np.savetxt('Y_mesh_'+str(ny)+'.txt',Yt)
elif p=='s':
np.savetxt('X_mesh_'+str(nx)+'.txt',Xt)
2np.savetxt('Y_mesh_'+str(ny)+'.txt',Yt)
return [Xt,Yt]
