import numpy as np
import matplotlib.pyplot as plt 
import os 

Nx = 301 #Number of spatial grid points
Nt = int(5e2) #Number of temporal grid points
hbar = 1.0
L = 1.0 #Simulation box length 
m = 1.0
k = np.pi/L
dx = L/(Nx - 1) 
dt = 1e-7

#Initial Setup
x = np.arange(0,L+dx,dx)
psi0 = np.sqrt(2/L)*np.sin((np.pi/L)*x)
psi0[0] = 0
psi0[-1] = 0
psi_mat = np.zeros([Nt,Nx],dtype='complex')
psi_mat[0][:] = psi0

#Define potential profile (Gaussian)
mu, sigma = 1.0/2.0, 1.0/20.0
V = -1e4*np.exp(-(x-mu)**2/(2.0*sigma**2))

for t in range(0,Nt-1):
    for i in range(1,Nx-1):
        psi_mat[t+1][i] = psi_mat[t][i] + (1j/2)*(dt/dx**2)*(psi_mat[t][i+1] - 2*psi_mat[t][i] + psi_mat[t][i-1]) - 1j*dt*V[i]*psi_mat[t][i]

    normal = np.sum(np.abs(psi_mat[t+1]**2))*dx
    #normalize 
    psi_mat[t+1] /= np.sqrt(normal) 

path_psi_mat = os.path.abspath("")
np.save('schrodinger_uwu/psi_mat.npy', psi_mat)
np.save('schrodinger_uwu/x_coords.npy', x)
np.save('schrodinger_uwu/potential.npy', V)











