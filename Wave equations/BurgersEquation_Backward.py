# -*- coding: utf-8 -*-
"""
Created on Sat May 13 16:30:08 2023

@author: jimli
"""

#%%

import numpy as np
import matplotlib.pyplot as plt

 # Parameters

N = 101
dx = 0.1
dt = dx * 0.1
iterations = 3001
t = 0.
σ = 2.
χμ = 5.
η = 0.2

x = np.linspace(0, 10, N)
U = np.zeros(N)

 # Initiate velocity

U = np.exp(-((x - χμ) / σ) ** 2)

 
#%%

 # Calculate velocity in Line

for j in range(iterations):
    
    U_n = np.copy(U)
    
    for i in range(1, N - 1):
        
        # Choose numerical derivative method
        
        if U_n[i] >= 0:
            dxU = (U_n[i] - U_n[i - 1]) / dx
        else:
            dxU = (U_n[i+1] - U_n[i]) / dx
        
        dxxU = (U_n[i + 1] - 2 * U_n[i] + U_n[i - 1]) / (dx**2)
        
        U[i] = U_n[i] + dt * (-U_n[i] * dxU + η * dxxU)
    
    # Harmonic Boundary Conditions
    
    U[N-1] = U[1]
    U[0] = U[N-2]   
    
    t = t + dt



# Plot animation 
    
    if j % 100 == 0:  # Frames
        plt.figure(figsize=(10, 6)) 
        plt.plot(x, U, 'r', label='Velocity') 
        plt.xlabel('x')
        plt.ylabel('U(x, t)')
        plt.title(f't = {j * dt:.2f}')
        plt.grid(True)  
        plt.legend()  
        plt.ylim(0,1.1)
plt.show()

