 # Libraries

import numpy as np
import matplotlib.pyplot as plt

 # Constants (S.I)

G = (6.67430e-11)**1 # gravitational constant 
M = (5.972e24)**1    # mass of the Earth
g = 9.81 # gravitational acceleration
m = 420000 # mass of the object in orbit

 # Initial conditions

x0 = 0.0    # initial x position 
x1_0 = 7778.0  # initial y position 
y0 = 6371000 + 408000 # initial x velocity 
y1_0 = 0.0   # initial y velocity

 # Time step/Max steps

dt = 1 # time step
tmax = 60*60*24*10 # simulation time
nsteps = int(tmax/dt) # max steps

 # Blank arrays to store values

E_t = np.zeros(nsteps)
E_x = np.zeros(nsteps)
E_x1 = np.zeros(nsteps)
E_y = np.zeros(nsteps)
E_y1 = np.zeros(nsteps)

RK4_t = np.zeros(nsteps)
RK4_x = np.zeros(nsteps)
RK4_x1 = np.zeros(nsteps)
RK4_y = np.zeros(nsteps)
RK4_y1 = np.zeros(nsteps)


 # Import initial conditions

RK4_t[0] = 0
RK4_x[0] = x0
RK4_x1[0] = x1_0
RK4_y[0] =y0
RK4_y1[0] = y1_0

E_t[0] = 0
E_x[0] = x0
E_x1[0] = x1_0
E_y[0] = y0
E_y1[0] = y1_0

 # Define functions

def fx1(RK4_x, RK4_y):
    return - G*M*RK4_x/(RK4_x**2 + RK4_y**2)**(3/2)

def fx(RK4_x1):
    return RK4_x1

def fy1(RK4_x, RK4_y):
    return - G*M*RK4_y/(RK4_x**2 + RK4_y**2)**(3/2)

def fy(RK4_y1):
    return RK4_y1

 # Application of Euler method

for i in range(nsteps-1):

   
    E_x1[i+1] = E_x1[i]  - G*M*E_x[i]/(E_x[i]**2 + E_y[i]**2)**(3/2)*dt
    E_x[i+1] = E_x[i] + E_x1[i]*dt
    
    E_y1[i+1] = E_y1[i]  - G*M*E_y[i]/(E_x[i]**2 + E_y[i]**2)**(3/2)*dt
    E_y[i+1] = E_y[i] + E_y1[i]*dt
    
    E_t[i+1] = E_t[i] + dt
    
 # Application of RK4 method

for i in range(nsteps-1):

        k1x1 = dt * fx1(RK4_x[i], RK4_y[i])
        k1x = dt * fx(RK4_x1[i])
        k1y1 = dt * fy1(RK4_x[i], RK4_y[i])
        k1y = dt * fy(RK4_y1[i])

        k2x1 = dt * fx1(RK4_x[i] + 0.5*k1x*dt, RK4_y[i] + 0.5*k1y*dt)
        k2x = dt * fx(RK4_x1[i] + 0.5*k1x1*dt)
        k2y1 = dt * fy1(RK4_x[i] + 0.5*k1x*dt, RK4_y[i] + 0.5*k1y*dt)
        k2y = dt * fy(RK4_y1[i] + 0.5*k1y1*dt)

        k3x1 = dt * fx1(RK4_x[i] + 0.5*k2x*dt, RK4_y[i] + 0.5*k2y*dt)
        k3x = dt * fx(RK4_x1[i] + 0.5*k2x1*dt)
        k3y1 = dt * fy1(RK4_x[i] + 0.5*k2x*dt, RK4_y[i] + 0.5*k2y*dt)
        k3y = dt * fy(RK4_y1[i] + 0.5*k2y1*dt)

        k4x1 = dt * fx1(RK4_x[i] + k3x*dt, RK4_y[i] + k3y*dt)
        k4x = dt * fx(RK4_x1[i] + k3x1*dt)
        k4y1 = dt * fy1(RK4_x[i] + k3x*dt, RK4_y[i] + k3y*dt)
        k4y = dt * fy(RK4_y1[i] + k3y1*dt)

        RK4_x1[i+1] = RK4_x1[i] + (1/6)*(k1x1 + 2*k2x1 + 2*k3x1 + k4x1)
        RK4_x[i+1] = RK4_x[i] + (1/6)*(k1x + 2*k2x + 2*k3x + k4x)
        RK4_y1[i+1] = RK4_y1[i] + (1/6)*(k1y1 + 2*k2y1 + 2*k3y1 + k4y1)
        RK4_y[i+1] = RK4_y[i] + (1/6)*(k1y + 2*k2y + 2*k3y + k4y)

        RK4_t[i+1] = RK4_t[i] + dt
        
 # Energy

E_T = 1/2 * m * (E_x1**2 + E_y1**2)
E_U = - G*M *m/ (np.sqrt(E_x**2 + E_y**2))
E_E = 1/2 * m * (E_x1**2 + E_y1**2) - G*M *m/ (np.sqrt(E_x**2 + E_y**2))


RK4_T = 1/2 * m * (RK4_x1**2 + RK4_y1**2)
RK4_U = - G*M *m/ (np.sqrt(RK4_x**2 + RK4_y**2))
RK4_E = 1/2 * m * (RK4_x1**2 + RK4_y1**2) - G*M *m/ (np.sqrt(RK4_x**2 + RK4_y**2))

 # Conver to Km

E_x /= 1000
E_y /= 1000
E_x1 /= 1000
E_y1 /= 1000       

RK4_x /= 1000
RK4_y /= 1000
RK4_x1 /= 1000
RK4_y1 /= 1000    

# Plot trajectory

plt.plot(E_x, E_y, linestyle='--',label='Euler',color='black')
plt.plot(RK4_x, RK4_y, linestyle='--',label='RK4',color='red')
plt.axhline(y=0, color='black', linestyle='-')
plt.axvline(x=0, color='black', linestyle='-')
plt.xlabel('x in kmeters')
plt.ylabel('y in kmeters')
plt.xlim(-np.max(np.abs(E_x)), np.max(np.abs(E_x)))
plt.ylim(-np.max(np.abs(E_y)), np.max(np.abs(E_y)))
plt.title('Orbit around the Earth')
plt.show()


# Plot Total mechanical energy

plt.plot(E_t, E_E, label='Total Energy in Euler',color='black')
plt.plot(E_t, E_U, label='Potential Energy in Euler',color='green')
plt.plot(E_t, E_T, label='Kinetic Energy in Euler',color='yellow')
plt.plot(RK4_t, RK4_E, label='Total Energy in RK4',color='red')
plt.plot(RK4_t, RK4_U, label='Potential Energy in RK4',color='blue')
plt.plot(RK4_t, RK4_T, label='Kinetic Energy in RK4',color='orange')
plt.legend(loc='upper left')
plt.xlabel('t in seconds')
plt.ylabel('E in joules')