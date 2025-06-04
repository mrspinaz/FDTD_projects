import matplotlib.pyplot as plt 
from matplotlib.animation import FuncAnimation, PillowWriter
import numpy as np
import plotly.express as px
from matplotlib.collections import LineCollection


psiR_mat = np.load('schrodinger_uwu/psiR_mat.npy')
psiI_mat = np.load('schrodinger_uwu/psiI_mat.npy')

prob_mat = np.load('schrodinger_uwu/prob_mat.npy')
V = np.load('schrodinger_uwu/potential.npy')
x = np.load('schrodinger_uwu/x_coords.npy')

V[V > -100] = np.nan
V *= 0.5
x /= np.max(x)

phase_mat = np.rad2deg(np.arctan2(psiI_mat, psiR_mat)) 
phase_mat[phase_mat < 0] += 360

fig=plt.figure()

axis = plt.axes(xlim =(0, 1), 
                ylim =(-10, 15))



line1, = axis.plot([], [], label='Real($\psi$)')
line3, = axis.plot([], [], label='Im($\psi$)')
line4, = axis.plot([], [], label='|$\psi^2$|')
line2, = axis.plot([], [], linewidth=3, label='_hidden')
line5 = LineCollection([], cmap=plt.cm.hsv,linewidth=4)
axis.add_collection(line5)


def init():
    line1.set_data([], [])
    line3.set_data([], [])
    line4.set_data([], [])
    line2.set_data([], [])
    line5.set_segments([])
    line5.set_array([])
    return line1, line2, line3, line4, line5,

def animate(i):
    y1 = psiR_mat[i*100][:]
    y2 = psiI_mat[i*100][:]
    y3 = prob_mat[i*100][:]
    

    line1.set_data(x,y1)
    line3.set_data(x,y2)
    line4.set_data(x,y3)
    line2.set_data(x,V*1e-3)

    #setup for phase color plotting 
    points = np.array([x,y3]).T.reshape(-1,1,2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    colors = phase_mat[i*100][:-1]
    line5.set_segments(segments)
    line5.set_array(colors)

    return line1, line2, line3, line4, line5,

anim = FuncAnimation(fig, animate, init_func=init, frames=2*999, interval=1, blit=True)
plt.xlabel('x/L',fontsize=15)
plt.show()

#print('Saving')
#writer = PillowWriter(fps=100)
#anim.save('schrodinger_uwu/test.gif', writer = writer)


