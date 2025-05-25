import matplotlib.pyplot as plt 
from matplotlib.animation import FuncAnimation, PillowWriter
import numpy as np

fig=plt.figure()

axis = plt.axes(xlim =(0, 1), 
                ylim =(-10, 20))

line1, = axis.plot([], [])
line2, = axis.plot([], [])
line3, = axis.plot([], [])
line4, = axis.plot([], [])



psiR_mat = np.load('schrodinger_uwu/psiR_mat.npy')
psiI_mat = np.load('schrodinger_uwu/psiI_mat.npy')
prob_mat = np.load('schrodinger_uwu/prob_mat.npy')
V = np.load('schrodinger_uwu/potential.npy')
x = np.load('schrodinger_uwu/x_coords.npy')

dx = 1.0/300
normal = sum(psiI_mat[5000][:]**2 + psiR_mat[5000][:]**2)*dx
print(normal)



def init():
    line1.set_data([], [])
    line2.set_data([], [])
    line3.set_data([], [])
    line4.set_data([], [])
    return line1, line2, line3, line4,

def animate(i):
    y1 = psiR_mat[i*100][:]
    y2 = psiI_mat[i*100][:]
    y3 = prob_mat[i*100][:]
    line1.set_data(x,y1)
    line2.set_data(x,V*1e-3)
    line3.set_data(x,y2)
    line4.set_data(x,y3)
    return line1, line2, line3, line4,

anim = FuncAnimation(fig, animate, init_func=init, frames=2*999, interval=1, blit=True)
plt.show()

#writer = PillowWriter(fps=60)
#anim.save('schrodinger_uwu/test.gif', writer = writer)


