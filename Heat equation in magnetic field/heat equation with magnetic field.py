# -*- coding: utf-8 -*-
"""
Created on Fri Apr 28 15:49:47 2023

@author: jimli
"""

import numpy as np
import matplotlib.pyplot as plt

 # Parameters
 
N = 11 # grid resolution
ke = 1.0 # thermal coefficient constant
Cv = 1.0
omega_b = 5.0
λ = 1.0
iterations = 4000 # number of iterations (time = iterations/10000 sec)
t=0.0

 # Create arrays to store values
 
x = np.linspace(0, 1, N)
y = np.linspace(0, 1, N)


T=np.zeros((N, N))
Q = np.zeros((N, N))


Bxx = np.zeros((N, N))
Byy = np.zeros((N, N))
Bxy = np.zeros((N, N))
Byx = np.zeros((N, N))
B = np.zeros((N, N))

Kxx = np.zeros((N, N))
Kyy = np.zeros((N, N))
Kxy = np.zeros((N, N))
Kyx = np.zeros((N, N))

 # spatial steps

x_min=0.0
x_max=1.0

y_min=0.0
y_max=1.0

dx=((x_max-x_min)/(N-1))
dy=((y_max-y_min)/(N-1))
dt = dx*dy*0.01 # time step


 # Calculate T using the Gauss-Seidel method
 
for k in range(0,iterations-1):
    
    for i in range(1, N-1):
        for j in range(1, N-1):
            
            xval = x[i]
            yval = y[j]
            
            B[i,j] = ((x[i]**3-x[i]**2)**2)*((1.1-2*y)[j]**2)+((3*x[i]**2-2*x[i])**2)*((1.1*y[j]-y[j]**2)**2)
            
            Byy[i][j] = ((3*x[i]**2-2*x[i])**2)*((1.1*y[j]-y[j]**2)**2)/ B[i,j]
            Bxx[i][j] = (((x[i]**3-x[i]**2)**2)*((1.1-2*y[j])**2))/ B[i,j]
            Byx[i][j] =  -(3*x[i]**2-2*x[i])*(1.1*y[j]-y[j]**2)*(x[i]**3-x[i]**2)*(1.1-2*y[j])/ B[i,j]
            
            
            dxBxx = (Bxx[i+1,j]-Bxx[i-1,j])/(2*dx)
            dyBxx = (Bxx[i,j+1]-Bxx[i,j-1])/(2*dy)
            dxByx = (Byx[i+1,j]-Byx[i-1,j])/(2*dx)
            dyByx = (Byx[i,j+1]-Byx[i,j-1])/(2*dy)
            dxB22 = (Byy[i+1,j]-Byy[i-1,j])/(2*dx)
            dyByy = (Byy[i,j+1]-Byy[i,j-1])/(2*dy)
            
            dxT = (T[i+1,j]-T[i-1,j])/(2*dx)
            dyT = (T[i,j+1]-T[i,j-1])/(2*dy)
            
            dxxT=(T[i+1,j]-2*T[i,j]+T[i-1,j])/dx**2
            dyyT=(T[i,j+1]-2*T[i,j]+T[i,j-1])/dy**2
            dxyT = (T[i+1,j+1]-T[i+1,j-1]-T[i-1,j+1]+T[i-1,j-1])/(4*dx*dy)
            
            
            Q[i][j] = ke * Cv * (1 + 25* Bxx[i][j])*dxxT \
                + ke *  Cv *  25*dxBxx*dxT \
                + ke * Cv *  25*Bxy[i][j]*dxyT \
                + ke * Cv *  25*dxByx*dyT \
                + ke * Cv *  25*dyByx*dxT \
                + ke * Cv *  25*Byx[i][j]*dxyT \
                + ke * Cv *  25*dyByy*dyT \
                + ke * Cv * (1 +  25*Byy[i][j])*dyyT \
                + ke * Cv * 1e3*np.exp(-((x[i]-0.5)**2 + (y[j]-0.5)**2)/(0.1)**2)*np.exp(-λ*t)

    for i in range(1,N-1):
        T[i,N-1] = T[i,N-2] - dy*T[i,N-2]**4

    for i in range (1, N-1):
        for j in range(1, N-1): 
            
                T[i,j]=T[i,j]+Q[i,j]*dt
         
    t=t+dt
      
 # Plot the temperature as a 2D surface
plt.style.use('dark_background')
plt.figure()
W, Z = np.meshgrid(np.arange(0, N), np.arange(0, N))
imageNbi = iterations-1
plt.contourf(Z/N, W/N, T, 50, cmap = plt.cm.inferno)
plt.colorbar()
plt.xlabel("x")
plt.ylabel("y")
plt.title("Temperature for = " + str(round(imageNbi*dt,2)) + " s")