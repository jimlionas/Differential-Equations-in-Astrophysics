# -*- coding: utf-8 -*-
"""
Created on Wed May 17 12:25:07 2023

@author: jimli
"""

#%%

import numpy as np
import matplotlib.pyplot as plt

 # Parameters

N = 1001
dx = 0.01
dt = dx * 0.1
iterations = 3001
t = 0.
γ = 5/3
η = 0.01

x = np.linspace(0, 10, N)
U = np.zeros(N)
P = np.ones(N)
ρ = np.ones(N)
e = np.zeros(N)
M = np.zeros(N)
cs = np.zeros(N)

 # Initiate velocity
 
U[0]=3.0
cs[0] =  np.sqrt(γ * P[0] /  ρ[0])
for i in range(0,N):
    e[i] = 1/2 * ρ[i] * U[i]**2 + (1/(γ-1))* P[i] 

 
#%%

 # Calculate velocity in Line

for j in range(iterations):
    
    U_n = np.copy(U)
    P_n = np.copy(P)
    ρ_n = np.copy(ρ)
    e_n = np.copy(e)
    M_n = np.copy(M)
    
    for i in range(1, N - 1):
        
        # Calculate derivatives
        
        dxxU = (U_n[i + 1] - 2 * U_n[i] + U_n[i - 1]) / (dx**2)
        
        dxU = (U_n[i+1] - U_n[i - 1]) /  (2*dx)
        
        dxU_b = (U_n[i] - U_n[i - 1]) /  (dx)
        
        dxρ = (ρ_n[i] - ρ_n[i - 1]) /  (dx)
        
        dxρ_c = (ρ_n[i+1] - ρ_n[i - 1]) /  (2*dx)
        
        dxe = (e_n[i] - e_n[i - 1]) /  (dx)
        
        dxP = (P_n[i+1] - P_n[i - 1]) /  (2*dx)
        
        # Calculate ρ,e,U,P 
        
        ρ[i] = ρ_n[i] + dt * (- ρ_n[i] * dxU - U_n[i] * dxρ)
        
        U[i] = U_n[i] + dt * (- U_n[i] * dxU_b - 1/ρ_n[i] * dxP + η * dxxU)
        
        #e[i] = e_n[i] + dt * (- 3/2 * ρ_n[i] * U_n[i]**2 * dxU - 1/2 * U_n[i]**3 * dxρ_c \
         #                     - γ/(γ-1) * U_n[i] * dxP  - γ/(γ-1) * P_n[i] * dxU)
        
        e[i] = e_n[i] + dt * (-e_n[i] * dxU - U_n[i] *dxe - P_n[i] * dxU - U_n[i] * dxP )
        
        P[i] = (e[i] - 1/2 * ρ[i] * U[i]**2) * (γ-1)
        
        M[i] = U[i] / np.sqrt(γ * P[i] / ρ[i])
        
        cs[i] = np.sqrt(γ * P[i] / ρ[i])
       
    
    t = t + dt

#%%

# Plot animation 
    
    if j % 1000 == 0:  # Frames
        plt.figure(figsize=(10, 6))  
        plt.plot(x, U, 'r', label='Velocity')  
        plt.xlabel('x')
        plt.ylabel('U(x, t)')
        plt.title(f't = {j * dt:.2f}')
        plt.grid(True) 
        plt.legend()  


    if j % 1000 == 0:  # Frames
        plt.figure(figsize=(10, 6))  
        plt.plot(x, P, 'yellow', label='Pressure')  
        plt.xlabel('x')
        plt.ylabel('P(x, t)')
        plt.title(f't = {j * dt:.2f}')
        plt.grid(True) 
        plt.legend()  

    if j % 1000 == 0:  # Frames
        plt.figure(figsize=(10, 6))  
        plt.plot(x, ρ , 'blue', label='Density')  
        plt.xlabel('x')
        plt.ylabel('ρ(x, t)')
        plt.title(f't = {j * dt:.2f}')
        plt.grid(True) 
        plt.legend() 

    if j % 1000 == 0:  # Frames
        plt.figure(figsize=(10, 6))  
        plt.plot(x, M, 'black', label='Mach number')  
        plt.xlabel('x')
        plt.ylabel('M(x, t)')
        plt.title(f't = {j * dt:.2f}')
        plt.grid(True) 
        plt.legend()  
        
        if j % 1000 == 0:  # Frames
            plt.figure(figsize=(10, 6))  
            plt.plot(x, cs, 'green', label='speed of sound')  
            plt.xlabel('x')
            plt.ylabel('Cs(x, t)')
            plt.title(f't = {j * dt:.2f}')
            plt.grid(True) 
            plt.legend()  
plt.show()

