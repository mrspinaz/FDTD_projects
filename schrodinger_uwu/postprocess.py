import matplotlib.pyplot as plt 
from matplotlib.animation import FuncAnimation, PillowWriter
import numpy as np

fig=plt.figure()

axis = plt.axes(xlim =(0, 1), 
                ylim =(-10, 20))

line1, = axis.plot([], [])
line2, = axis.plot([], [])
line3, = axis.plot([], [])



psiR_mat = np.load('schrodinger_uwu/psiR_mat.npy')
psiI_mat = np.load('schrodinger_uwu/psiI_mat.npy')
V = np.load('schrodinger_uwu/potential.npy')
x = np.load('schrodinger_uwu/x_coords.npy')

def init():
    line1.set_data([], [])
    line2.set_data([], [])
    line3.set_data([], [])
    return line1, line2, line3,

def animate(i):
    y1 = psiR_mat[i*10][:]
    y2 = psiI_mat[i*10][:]

    line1.set_data(x,y1)
    line2.set_data(x,V*1e-4)
    line3.set_data(x,y2)
    return line1, line2, line3,

anim = FuncAnimation(fig, animate, init_func=init, frames=1000, interval=1, blit=True)
plt.show()
writer = PillowWriter(fps=30)
anim.save('schrodinger_uwu/test.gif', writer = writer)


