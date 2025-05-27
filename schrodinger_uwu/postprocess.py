import matplotlib.pyplot as plt 
from matplotlib.animation import FuncAnimation, PillowWriter
import numpy as np

fig=plt.figure()

axis = plt.axes(xlim =(0, 1), 
                ylim =(-10, 15))



line1, = axis.plot([], [], label='Real($\psi$)')
line3, = axis.plot([], [], label='Im($\psi$)')
line4, = axis.plot([], [], label='|$\psi^2$|')
line2, = axis.plot([], [],label='_hidden')



psiR_mat = np.load('schrodinger_uwu/psiR_mat.npy')
psiI_mat = np.load('schrodinger_uwu/psiI_mat.npy')
prob_mat = np.load('schrodinger_uwu/prob_mat.npy')
V = np.load('schrodinger_uwu/potential.npy')
x = np.load('schrodinger_uwu/x_coords.npy')
V[V > -1000] = np.nan
x /= np.max(x)



def init():
    line1.set_data([], [])
    line3.set_data([], [])
    line4.set_data([], [])
    line2.set_data([], [])
    return line1, line2, line3, line4,

def animate(i):
    y1 = psiR_mat[i*100][:]
    y2 = psiI_mat[i*100][:]
    y3 = prob_mat[i*100][:]
    line1.set_data(x,y1)
   
    line3.set_data(x,y2)
    line4.set_data(x,y3)
    line2.set_data(x,V*1e-3)
    return line1, line2, line3, line4,

anim = FuncAnimation(fig, animate, init_func=init, frames=3*999, interval=1, blit=True)
plt.xlabel('x/L',fontsize=15)
plt.legend(fontsize=15)
plt.show()

#writer = PillowWriter(fps=60)
#anim.save('schrodinger_uwu/test.gif', writer = writer)


