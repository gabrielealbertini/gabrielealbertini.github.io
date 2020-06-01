#!/usr/bin/env python2
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.gridspec as gridspec
#import cv2

from postprocess.at_time import *

from postprocess.space_time import *

from utilities.latexify import *

latexify_jmps(plt,fig_width=2.75,fig_height=2.85)
#latexify(plt,style='sans',fig_width=4,fig_height=2.25)

plt.rcParams['animation.html'] = 'html5'

plt.rcParams['xtick.labelsize']= 5
plt.rcParams['ytick.labelsize']= 5

bname='stochastic_50'
group='interface'
#plot_space_time(bname,group,FieldId('cohesion',0))
#plt.show()

N=40;
tend = 0.00205
Tcr=0.00187
taupm=1e6
L=0.1

def data_gen():
    t = data_gen.t
    cnt = 0
    dt= tend/N
    while t < tend+dt:
        cnt+=1
        t += dt
        if t>0.9*Tcr:
            dt=tend/N/15
        if t>Tcr:
            dt=tend/N/30
        print(cnt,t)
        x,tau = get_at_time(bname, group, FieldId('cohesion',0),time=t)
        x,vel = get_at_time(bname, group, FieldId('top_velo',0),time=t)
        vel +=1e-8
        # adapted the data generator to yield 
        yield t, x, tau, vel

data_gen.t = 0

# create a figure with three subplots 
fig = plt.figure()
spec0 = gridspec.GridSpec(ncols=1,nrows=1,left=0.02,right=0.5,top=1,bottom=0.68,hspace=0.2)
spec1 = gridspec.GridSpec(ncols=2,nrows=1,left=0.15,right=0.96,top=0.97,bottom=0.8,hspace=0.2)
spec2 = gridspec.GridSpec(ncols=1,nrows=2,left=0.15, right=0.96,top=0.65,bottom=0.15,hspace=0.2)
ax0 = fig.add_subplot(spec0[0, 0])
ax1 = fig.add_subplot(spec1[0, 1])#, sharex=ax)
ax2 = fig.add_subplot(spec2[0, 0])
ax3 = fig.add_subplot(spec2[1, 0])

ax0.axis('off')
img = plt.imread('fig01_model_inset_v0.png')
print(img.shape)
im = ax0.imshow(img)
#define labels and limits
ax1.set_xlabel(r'$\mathrm{time}$ $t/T$')
ax1.set_ylabel(r' $\tau_0/\langle\tau_p\rangle$')#$\mathrm{applied \ stress}$
ax2.set_ylabel(r'$\mathrm{stress}$ $\tau/\langle\tau_p\rangle$')
ax3.set_xlabel(r'$\mathrm{position}$ $x/L$')
ax3.set_ylabel(r'$\mathrm{slip \ rate}$ $\partial u/\partial t$')

ax1.set_xticks(np.arange(0,2,0.5))
ax1.set_xlim([0,1.25])
ax1.set_ylim([0,1])

ax2.set_xlim([0,1])
ax2.set_yticks(np.arange(0,2.5,0.5))
ax2.set_ylim([0.0,2.0])
ax2.set_xticklabels([])
ax2.set_xticks(np.arange(0,1.1,0.25))

ax3.set_xticks(np.arange(0,1.1,0.25))
ax3.set_xlim([0,1])
ax3.set_ylim([1e-6,1e1])
ax3.set_yscale('log')
ax3.get_yaxis().set_tick_params(which='minor', size=0.75)


# intialize two line objects (one in each axes)
x00,tau00 = get_at_time(bname, group, FieldId('cohesion',0),time=0)

with open('data/stochastic_50-datamanager-files/stochastic_50.x.random') as f:
    x = f.read()
    x = np.array([float(xi) for xi in x.split(',')])
with open('data/stochastic_50-datamanager-files/stochastic_50.tauc.random') as f:
    taup = f.read()
    taup = np.array([float(xi) for xi in taup.split(',')])

ax2.plot(x/L,taup/taupm,color=np.array([1,1,1])*0.7)

ax2.annotate(r'$\mathrm{(Albertini \ et \  al. \ 2020)}$', xy=(-0.175, -3.55), annotation_clip=False,)#, xycoords='data'
line0, = ax1.plot([], [], color=matlab_colors[6])
line1, = ax2.plot(x00/L, tau00/taupm, color=matlab_colors[0])
line2, = ax3.plot([], [], color=matlab_colors[4])


plt.savefig('fig03_movie_web.png',dpi=300)
plt.savefig('fig03_movie_web.pdf',dpi=600)


line = [line0, line1, line2]

# initialize the data arrays 
time, space, glstress, stress, vel = [],[],[],[],[];
time.append(0.0)
glstress.append(np.min(taup))
def run(data):
    # update the data
    t, space, stress, velo = data
    time.append(t)
    glstress.append(np.mean(stress))

    # update the data of  line objects
    line[0].set_data(np.array(time)/Tcr,np.array(glstress)/taupm)
    line[1].set_data(space/L,stress/taupm)
    line[2].set_data(space/L,velo)

    return line

    
ani = animation.FuncAnimation(fig, run, data_gen, blit=True, interval=100,
                              repeat=False,save_count=sys.maxsize)


ani.save('fig03_movie_web.gif',writer='imagemagick',dpi=300)

#ani.save('fig03_movie.mp4',codec="libx264",bitrate=-1,dpi=300)
#plt.show()
