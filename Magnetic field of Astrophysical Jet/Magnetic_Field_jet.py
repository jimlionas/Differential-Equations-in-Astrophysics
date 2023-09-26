

import numpy as np
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import matplotlib

 # Values for B0,B1,R0
 
B0 = 1 
B1 = 20
R0 = 1


 # Initialize grid
 
R_min = 0.0
R_max = 5.0
Z_min = 0.0
Z_max = 5.0

N_R = 101 #R resolution
N_Z = 101 #Z resolution

dR = (R_max - R_min) / (N_R - 1)
dZ = (Z_max - Z_min) / (N_Z - 1)


 # Create arrays to store values

I = np.zeros((N_R, N_Z))
Psi = np.zeros((N_R, N_Z))
Result=np.zeros((3,N_R*N_Z))
R=np.zeros(N_R)
Z=np.zeros(N_Z)

 # Setting the values of R and Z

for i in range(0,N_R):
    R[i] = i * dR
    
for i in range(0,N_Z):
    Z[i] = i * dZ


 # Set the boundary conditions for I and Ψ on the grid


Psi[:,0] = 0 # left
Psi[N_R-1,:] = Psi[N_R-2,i]  # up
Psi[:,N_R-1] =  B0 * R_max**2 / 2  # right

for j in range (0,N_R):
    Psi[0,j]= B0 * R[j]**2 / 2 # down
    

I = B1 *  R0 * np.sin(np.pi * 2 * Psi / (B0 *R_max**2))
dI = B1 * R0 * ((2 * np.pi) / (B0 * R_max**2)) * np.cos((2*np.pi*Psi)/(B0*R_max**2))

 # number of iterations
iterations = 5000

L=np.zeros((N_R,N_Z))
for k in range(0, iterations):
    
    #In this loop we evaluate the derivatives
    for i in range (1, N_Z-1):
        for j in range(1, N_R-1):
            dRRPsi=(Psi[i,j+1]-2*Psi[i,j]+Psi[i,j-1])/dR**2
            dZZPsi=(Psi[i+1,j]-2*Psi[i,j]+Psi[i-1,j])/dZ**2
            dRPsi = (Psi[i+1,j] - Psi[i,j]) / dR
        
            I[i,j] = B1 *  R0 * np.sin(np.pi * 2 * Psi[i,j] / (B0 *R_max**2))
            dI[i,j] = B1 * R0 * ((2 * np.pi) / (B0 * R_max**2)) * np.cos((2*np.pi*Psi[i,j])/(B0*R_max**2))
            
            L[i,j]=dRRPsi+dZZPsi +(1/R[i])*dRPsi + I[i,j]*dI[i,j]
            
         
           
    
    #In this loop we evaluate the new value of Psi
            
    for i in range (1, N_R-1):
        for j in range(1, N_Z-1): 
            
            Psi[i,j]=Psi[i,j]+L[i,j]*dR*dZ*0.1
            
 # In this loop we save the results in a single array
            
    for  j in range(1,N_R):
        Psi[N_Z-1,j] = Psi[N_Z-2,j]
for i in range(0,N_R):
    for j in range(0, N_Z):
        k=i*N_R+j
        Result[0, k]=R[i]
        Result[1, k]=Z[j]
        Result[2, k]=Psi[j,i]

X= Result[0, :].reshape(N_R,N_Z)
Y= Result[1, :].reshape(N_R,N_Z)
Z= Result[2, :].reshape(N_R,N_Z)

 # plot the magnetic field

fig = plt.figure(dpi=1000)
ax1 = fig.add_subplot(1,1,1, aspect=1, xlim=[R_min, R_max], ylim=[Z_min, Z_max])
ax1.set_xlabel('R')
ax1.set_ylabel('Z')
ax1.set_title('Magnetic Field')
cf = ax1.contourf(X,Y,Z)
cbar = fig.colorbar(cf, ax=ax1)
cbar.set_label('Magnetic field strength')
plt.savefig('magnetic_field.png')


plt.savefig('Laplacian.png')

 # Plot the magnetic field lines
 
R = np.linspace(R_min, R_max, N_R)
Z = np.linspace(Z_min, Z_max, N_Z)
R, Z = np.meshgrid(R, Z)


fig, ax = plt.subplots()
ax.set_xlabel('R')
ax.set_ylabel('Z')
ax.set_title('Magnetic Field Lines')
contour = ax.contour(R, Z, Psi, 100, cmap='Reds')
cbar = plt.colorbar(contour)
cbar.set_label('Magnetic field lines')
plt.show()

 # plot the current I
 
I[i,j] = B1 *  R0 * np.sin(np.pi * 2 * Psi[i,j] / (B0 *R_max**2))
 
fig, ax = plt.subplots()
ax.set_xlabel('R')
ax.set_ylabel('Z')
ax.set_title('Current')
contourf = ax.contourf(R, Z, I, 100, cmap='PuBu')
cbar = plt.colorbar(contourf)
cbar.set_label('Current intesity')
plt.show()
