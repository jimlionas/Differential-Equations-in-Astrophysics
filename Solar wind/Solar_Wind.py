# -*- coding: utf-8 -*-
"""
Created on Mon May 22 10:52:11 2023

@author: jimli
"""

import numpy as np
import matplotlib.pyplot as plt

 # Parameters

N = 101
dr = 0.1
dt = dr * 0.01
iterations = 100001
t = 0.
K = 1.
Μ = 5.
G = 1.
η = 0.1
r_min = 1.
r_max = 11


r = np.linspace(r_min, r_max, N)
U = np.zeros(N)
P = np.ones(N) * 0.01
ρ = np.ones(N) * 0.01
cs = np.zeros(N)

 # Initiate velocity
 
ρ[0]=1.
P[0]=1.
U[0] = 0.19
cs[0] = np.sqrt(1 / ρ[0])

 
#%%

 # Calculate velocity in Line

for j in range(iterations):
    
    U_n = np.copy(U)
    P_n = np.copy(P)
    ρ_n = np.copy(ρ)
    
    for i in range(1, N - 1):
        
        # Calculate derivatives
        
        drrU = (U_n[i + 1] - 2 * U_n[i] + U_n[i - 1]) / (dr**2)
        
        drU = (U_n[i ] - U_n[i - 1]) /  (dr)
        drU_c = (U_n[i + 1] - U_n[i - 1]) /  (2*dr)
        
        drρ = (ρ_n[i] - ρ_n[i - 1]) /  (dr)
            
        # Calculate ρ,e,U,P 
        
        ρ[i] = ρ_n[i] + dt * ((- 2/r[i]) * ρ_n[i] * U_n[i]  - drρ * U_n[i] - drU_c * ρ_n[i])
        
        U[i] = U_n[i] + dt * (-U_n[i] * drU - (K/ρ_n[i])* drρ - (Μ * G) / r[i]**2 + η * ((2/r[i])*drU+drrU))
        
        cs[i] = np.sqrt(1 / ρ[i])

        
        
        
        P[i] = K * ρ[i]
        
    
    t = t + dt
    U[N-1] = U[N-2]
    ρ[N-1] = ρ[N-2]
    P[N-1] = P[N-2]

#%%


# Plot animation 
    
    if j % 1000 == 0:  # Frames
        plt.figure(figsize=(10, 6))
        plt.plot(r, U, 'r', label='Velocity')
        plt.xlabel('r')
        plt.ylabel('U(x, t)')
        plt.title(f't = {j * dt:.2f}')
        plt.grid(True)
        plt.axvline(x=2.5, color='black', linestyle='--', label='r = 2.5')
        
plt.legend()
    