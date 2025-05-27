import numpy as np
import matplotlib.pyplot as plt 

def gaussian_wavepacket(x, center, sigma, p0):

    gaussian = np.exp(-(x-center)**2/(2*sigma**2))

    return gaussian

def fermi(x,a,b,c,d):
    fermi_func = (1.0*a)/(1.0*b + np.exp((x + c)/d))
    return fermi_func

def u_eye(x,width,depth,position):

    f1 = fermi(x, depth, 0.1,position, 0.05)
    f2 = fermi(x, depth, 0.1,-width + position, 0.05)

    u_eye = f1-f2
    return f1 - f2


Nx = 401 #Number of spatial grid points
Nt = int(3e5) #Number of temporal grid points
hbar = 1.0
L = 3.0 #Simulation box length 
m = 1.0
k = np.pi/L
dx = L/(Nx - 1) 
dt = 1e-7
c1 = (hbar*dt)/(2*m*dx**2)
c2 = dt/hbar

#Initial Setup
x = np.arange(0,L+dx,dx)

psiR_mat = np.zeros([Nt,Nx])
psiI_mat = np.zeros([Nt,Nx])
prob_mat = np.zeros([Nt,Nx])

#Set up initial conditions:
#psi0 = np.sqrt(2/L)*np.sin((np.pi/L)*x)
gw = gaussian_wavepacket(x, center=2.3, sigma=0.05, p0=-80.0) #Get gaussian envelope

psiR0 = np.cos(-80*x)*gw #Gaussian envelope 
psiI0 = np.sin(-80*x)*gw #Gaussian envelope
normal_init = sum(psiR0**2 + psiI0**2)*dx #normalize

psiR0 /= np.sqrt(normal_init) #normalize
psiI0 /= np.sqrt(normal_init) #normalize

psiR0[0] = 0 #BC
psiI0[0] = 0 #BC
psiR0[-1] = 0 #BC
psiI0[-1] = 0 #BC

psiR_mat[0][:] = psiR0
psiI_mat[0][:] = psiI0
prob_mat[0][:] = psiR_mat[0][:]**2 + psiI_mat[0][:]**2
#psiI initially should be zeros for infinite well ground state.

#Define potential profile
mu, sigma = 1.0/2.0, 1.0/20.0
V = -1e4*np.exp(-(x-mu)**2/(2.0*sigma**2))

#eyes
left_eye = u_eye(x*1.5, width=0.5, depth=300.0, position=-0.4)
right_eye = u_eye(x*1.5, width=0.5, depth=300.0, position=-2.4)

#mouth
left_part = u_eye(x*1.5, width=0.25, depth=300.0, position=-1.3)
right_part = u_eye(x*1.5, width=0.25, depth=300.0, position=-1.7)
mouth = right_part + left_part

V = -1*(left_eye + right_eye + mouth + 800)
V_scaled = V*1e-3
V_scaled[V_scaled > 0] = np.nan
plt.plot(x,-1*V_scaled)
plt.xlim(0,3)
plt.ylim(-15,20)
plt.show()



for t in range(0,Nt-1):
    
    psiI_mat[t+1][1:Nx-1] =  c1*(psiR_mat[t][2:Nx] - 2*psiR_mat[t][1:Nx-1] + psiR_mat[t][0:Nx-2]) - c2*V[1:Nx-1]*psiR_mat[t][1:Nx-1] + psiI_mat[t][1:Nx-1]
      
    psiR_mat[t+1][1:Nx-1] =  -c1*(psiI_mat[t+1][2:Nx] - 2*psiI_mat[t+1][1:Nx-1] + psiI_mat[t+1][0:Nx-2]) + c2*V[1:Nx-1]*psiI_mat[t+1][1:Nx-1] + psiR_mat[t][1:Nx-1]
    
    normal = np.sum(psiR_mat[t+1][1:Nx-1]**2 + psiI_mat[t+1][1:Nx-1]**2)*dx
    #normalize 
    psiR_mat[t+1][1:Nx-1] /= np.sqrt(normal)
    psiI_mat[t+1][1:Nx-1] /= np.sqrt(normal) 
    prob_mat[t+1][1:Nx-1] = (psiI_mat[t+1][1:Nx-1]**2 + psiR_mat[t+1][1:Nx-1]**2)/normal



np.save('schrodinger_uwu/psiR_mat.npy', psiR_mat)
np.save('schrodinger_uwu/psiI_mat.npy', psiI_mat)
np.save('schrodinger_uwu/prob_mat.npy', prob_mat)
np.save('schrodinger_uwu/x_coords.npy', x)
np.save('schrodinger_uwu/potential.npy', V)

print("Calculation Completed")












