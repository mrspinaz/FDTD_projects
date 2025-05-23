import matplotlib.pyplot as plt 
from matplotlib.animation import FuncAnimation, PillowWriter
import numpy as np

fig=plt.figure()

axis = plt.axes(xlim =(0, 1), 
                ylim =(-10, 20))

line1, = axis.plot([], [])
line2, = axis.plot([], [])


V = np.load('schrodinger_uwu/potential.npy')
psi_mat = np.load('schrodinger_uwu/psi_mat.npy')
x = np.load('schrodinger_uwu/x_coords.npy')

def init():
    line1.set_data([], [])
    line2.set_data([], [])
    return line1, line2,

def animate(i):
    y = np.abs(psi_mat[i*10][:])**2
    line1.set_data(x,y)
    line2.set_data(x,V*1e-4)
    return line1, line2,

anim = FuncAnimation(fig, animate, init_func=init, frames=5000, interval=1, blit=True)
plt.show()
writer = PillowWriter(fps=30)
anim.save('schrodinger_uwu/test.gif', writer = writer)


