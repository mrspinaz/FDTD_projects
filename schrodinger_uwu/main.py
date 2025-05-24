import numpy as np
import matplotlib.pyplot as plt 

Nx = 301 #Number of spatial grid points
Nt = int(1e4) #Number of temporal grid points
hbar = 1.0
L = 1.0 #Simulation box length 
m = 1.0
k = np.pi/L
dx = L/(Nx - 1) 
dt = 1e-7
c1 = (hbar*dt)/(2*m*dx**2)
c2 = dt/hbar

#Initial Setup
x = np.arange(0,L+dx,dx)
psi0 = np.sqrt(2/L)*np.sin((np.pi/L)*x)
psi0[0] = 0
psi0[-1] = 0
psiR_mat = np.zeros([Nt,Nx],dtype='float')
psiI_mat = np.zeros([Nt,Nx],dtype='float')
psiR_mat[0][:] = psi0
#psiI initially should be zeros for infinite well ground state.

#Define potential profile (Gaussian)
mu, sigma = 1.0/2.0, 1.0/20.0
V = x*0 #-1e4*np.exp(-(x-mu)**2/(2.0*sigma**2))

for t in range(0,Nt-1):
    for i in range(1,Nx-1):
        psiI_mat[t+1][i] =  c1*(psiR_mat[t][i+1] - 2*psiR_mat[t][i] + psiR_mat[t][i-1]) - c2*V[i]*psiR_mat[t][i] + psiI_mat[t][i]
        psiR_mat[t+1][i] =  -c1*(psiI_mat[t+1][i+1] - 2*psiI_mat[t+1][i] + psiI_mat[t+1][i-1]) + c2*V[i]*psiI_mat[t+1][i] + psiR_mat[t][i]

    normal = np.sum(psiR_mat[t+1][i]**2 + psiI_mat[t+1][i]**2)*dx
    #normalize 
    psiR_mat[i+1] /= np.sqrt(normal)
    psiI_mat[t+1] /= np.sqrt(normal) 

np.save('schrodinger_uwu/psiR_mat.npy', psiR_mat)
np.save('schrodinger_uwu/psiI_mat.npy', psiI_mat)
np.save('schrodinger_uwu/x_coords.npy', x)
np.save('schrodinger_uwu/potential.npy', V)











