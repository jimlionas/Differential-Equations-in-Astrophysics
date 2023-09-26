 # Libraries

import numpy as np
import matplotlib.pyplot as plt

 # Constants (S.I)

G = (6.67430e-11)**1 # gravitational constant 
M = (5.972e24)**1    # mass of the Earth
g = 9.81 # gravitational acceleration
m = 420000 # mass of the object in orbit

 # Initial conditions

r1_0 =  0.0      # initial x position 
r_0 = 6371000.0 + 408000.0  # initial y position 
omega_0 = 0.001135 # initial x velocity 
theta_0 = 0.0   # initial y velocity

 # Time step/Max steps

dt = 1 # time step
tmax = 60*60*24*10 # simulation time
nsteps = int(tmax/dt) # max steps

 # Blank arrays to store values

E_t = np.zeros(nsteps)
E_r = np.zeros(nsteps)
E_r1 = np.zeros(nsteps)
E_theta = np.zeros(nsteps)
E_omega = np.zeros(nsteps)

RK4_t = np.zeros(nsteps)
RK4_r = np.zeros(nsteps)
RK4_r1 = np.zeros(nsteps)
RK4_theta = np.zeros(nsteps)
RK4_omega = np.zeros(nsteps)


 # Import initial conditions

RK4_t[0] = 0
RK4_r[0] = r_0
RK4_r1[0] = r1_0
RK4_omega[0] = omega_0
RK4_theta[0] = theta_0

E_t[0] = 0
E_r[0] = r_0
E_r1[0] = r1_0
E_omega[0] = omega_0
E_theta[0] = theta_0

 # Define functions

def fr1(RK4_r,RK4_omega):
    return RK4_r * RK4_omega**2 - G*M/(RK4_r**2) 

def fr(RK4_r1):
    return RK4_r1

def fomega(RK4_r, RK4_omega, RK4_r1):
    return -2 * RK4_r1 * RK4_omega / RK4_r

def ftheta(RK4_omega):
    return RK4_omega

 # Application of Euler method

for i in range(nsteps-1):

   

   E_r1[i+1] = E_r1[i] + (E_r[i] * E_omega[i]**2 - G*M/((E_r[i]**2) )) * dt
   E_r[i+1] = E_r[i] + E_r1[i]*dt
    
   E_omega[i+1] = E_omega[i]  - (2 * E_r1[i] * E_omega[i]/E_r[i] ) * dt
   E_theta[i+1] = E_theta[i] + E_omega[i]*dt
    
   E_t[i+1] = E_t[i] + dt

 # Application of RK4 method

for i in range(nsteps-1):

    k1r1 =  dt * fr1(RK4_r[i], RK4_omega[i])
    k1r = dt * fr(RK4_r1[i])
    k1omega =  dt * fomega(RK4_r[i], RK4_omega[i], RK4_r1[i])
    k1theta = dt * ftheta(RK4_omega[i])

    k2r1 = dt *  fr1(RK4_r[i] + 0.5*k1r*dt, RK4_omega[i] + 0.5*k1omega*dt)
    k2r = dt *  fr(RK4_r1[i] + 0.5*k1r1*dt)
    k2omega = dt *  fomega(RK4_r[i] + 0.5*k1r*dt, RK4_omega[i] + 0.5*k1omega*dt, RK4_r1[i] + 0.5*k1r1*dt)
    k2theta = dt *  ftheta(RK4_omega[i] + 0.5*k1omega*dt)

    k3r1 = dt *  fr1(RK4_r[i] + 0.5*k2r*dt, RK4_omega[i] + 0.5*k2omega*dt)
    k3r = dt *  fr(RK4_r1[i] + 0.5*k2r1*dt)
    k3omega = dt *  fomega(RK4_r[i] + 0.5*k2r*dt, RK4_omega[i] + 0.5*k2omega*dt, RK4_r1[i] + 0.5*k2r1*dt)
    k3theta = dt *  ftheta(RK4_omega[i] + 0.5*k2omega*dt)

    k4r1 = dt *  fr1(RK4_r[i] + k3r*dt, RK4_omega[i] + k3omega*dt)
    k4r = dt *  fr(RK4_r1[i] + k3r1*dt)
    k4omega = dt *  fomega(RK4_r[i] + k3r*dt, RK4_omega[i] + k3omega*dt, RK4_r1[i] + k3r1*dt)
    k4theta = dt *  ftheta(RK4_omega[i] + k3omega*dt)

    RK4_r1[i+1] = RK4_r1[i] + (1/6)*(k1r1 + 2*k2r1 + 2*k3r1 + k4r1)
    RK4_r[i+1] = RK4_r[i] + (1/6)*(k1r + 2*k2r + 2*k3r + k4r)
    RK4_omega[i+1] = RK4_omega[i] + (1/6)*(k1omega + 2*k2omega + 2*k3omega + k4omega)
    RK4_theta[i+1] = RK4_theta[i] + (1/6)*(k1theta + 2*k2theta + 2*k3theta + k4theta)

    RK4_t[i+1] = RK4_t[i] + dt
    
 # Energy

E_E =(1/2)* m*E_r1**2 +(1/2)* m * E_r**2 * E_omega**2  - m *G*M / E_r
E_U = - G*M *m / E_r
E_T = (1/2)* m*E_r1**2 + (1/2)*E_r**2 * E_omega**2 

RK4_E =(1/2)* m*RK4_r1**2 +(1/2)* m * RK4_r**2 * RK4_omega**2  - m *G*M / RK4_r
RK4_U = - G*M *m / RK4_r
RK4_T = (1/2)* m*RK4_r1**2 + (1/2)*RK4_r**2 * RK4_omega**2 

 # Conver to Km

RK4_r /= 1000
E_r /= 1000

 # Plot trajectory


plt.axes(projection='polar')
plt.polar(E_theta,E_r,linestyle='--',color='black')
plt.polar(RK4_theta,RK4_r,linestyle='--',color='red')
plt.xlabel('Radius in (km)')
plt.title('Orbit around the Earth')
plt.show()

 # Plot r(t)

plt.plot(E_t, E_r, label='radius in Euler',color='black')
plt.plot(RK4_t, RK4_r, label='radius in RK4',color='red')
plt.xlabel('t in seconds')
plt.ylabel('radius in kmeters')
plt.legend(loc='upper left')
plt.show()

 # Plot Total mechanical energy 

plt.plot(E_t, E_E, label='Total Energy in Euler',color='black')
plt.plot(E_t, E_U, label='Potential Energy in Euler',linestyle='--',color='yellow')
plt.plot(E_t, E_T, label='Kinetic Energy in Euler',linestyle='--',color='green')
plt.plot(RK4_t, RK4_E, label='Total Energy in RK4',color='red')
plt.plot(RK4_t, RK4_U, label='Potential Energy in RK4',linestyle='--',color='blue')
plt.plot(RK4_t, RK4_T, label='Kinetic Energy in RK4',linestyle='--',color='orange')
plt.legend(loc='upper left')
plt.xlabel('t in seconds')
plt.ylabel('E in joules')