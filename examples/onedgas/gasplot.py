#!/usr/bin/python

"""
make an animation of the centre of mass positions, coloured by some physically relevant quantity.
can save movie file.

more details:
animate-phases.py -h
"""
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.collections as mc
import matplotlib.colors as colors
import matplotlib.cm as cmx
from matplotlib import animation
import sys
import os
from scipy.spatial import Delaunay as Del
import argparse
import signal
import polydata

dotsize=150

def weimueller(x,y,v1,v2,p=1,q=2,T=1):
    return v1*np.cos(q*x)-v2*np.cos(p*x-2*np.pi*y/T)

def parse_fname(fname):
    def get_digit(token):
        return float(fname.split(token)[1].split('_')[0])

    Lx, Ly, N, nm, v1, v2 = [get_digit(token) for token in ['Lx','Ly','N','np','v1','v2']]
    return Lx, Ly, N, nm, v1, v2

def vis2d(x,y,v1=0.,v2=0.,repx=10,repy=5,p=1,q=2,lw=2,colormap='hsv',npts=100):
    fig = plt.figure(figsize=(8,8))
    ax = fig.add_subplot(1, 1, 1, xticks=[], yticks=[])
    fig.subplots_adjust(0.05,0.05,0.95,0.95)
    fig.figurePatch.set_edgecolor('white')
    fig.figurePatch.set_facecolor('white')
#    ax.set_aspect('equal')
    ax.set_xlim([-Lx/2.,Lx/2.])
    ax.set_ylim([-Ly/2.,Ly/2.])
    ax.set_axis_off()

    if v1**2 + v2**2 > 0:
        onestepfactor = 2*np.pi/(Lx*1./repx)
        Tt = Ly*1./repy
        x,y = np.mgrid[-Lx/2.:Lx/2.:npts*1j, -Ly/2.:Ly/2.:npts*1j] 
        pot = weimueller(x,y,v1,v2,p=p*onestepfactor,q=q*onestepfactor,T=Tt)
        ax.pcolormesh(x,y,pot,cmap='RdBu')
    
    plotpts = plt.scatter(x,y,edgecolors='none',s=dotsize)
    return fig, plotpts

#############################################



# capture ctrl-c
def signal_handler(signal, frame):
    sys.exit(0)

signal.signal(signal.SIGINT, signal_handler)


## Parse arguments
parser = argparse.ArgumentParser(description='Real-time display of pt data.')

parser.add_argument('directory',
                    help='data folder to analyze; expected to have \'pts\' and \'raw\' subfolders')
parser.add_argument('--init',type=int,default=0,help='first frame')
parser.add_argument('--fin',type=int,default=-1,help='last frame')
parser.add_argument('--step',type=int,default=1,help='steps between frames')
parser.add_argument('--interval',type=int,default=5,help='pause between frames')
parser.add_argument('--phase',type=int,nargs='?',const=1,help='colour by angular dof, multiplied by optional argument')
parser.add_argument('--spin',type=float,default=0.,help='spin velocity to subtract from evolution')
parser.add_argument('--dotsize',type=float,default=150.,help='scatter plot dot size')
parser.add_argument('--persist',type=int,nargs='?',const=1,help='do not erase first frame')
parser.add_argument('-o', '--moviefile',help='if provided, saves output to movie filename')

args = parser.parse_args()

pd = polydata.PolyData(args.directory)

Lx, N, v1, v2 = [pd.params[i] for i in ["Lx","N","v1","v2"]]
try:
    repx, repy, p, q = [pd.params[i] for i in ["repx","repy","wmp","wmq"]]
except KeyError:
    repx = 10
    repy = 5
    p = 1
    q = 2

ptfiletemplate = args.directory+'/pts/sim%08d.pts'

data = []

for frame in np.arange(pd.lastframe)[::1000]:
    ptfile = ptfiletemplate % (frame)
    coords = np.loadtxt(ptfile)
    data.append(coords)

data = np.array(data)
print data.shape


# fig,plotpts = vis2d(x,y,v1,v2,repx,repy,p,q)
plt.figure()

factor = 2*np.pi/Lx
datauw = np.unwrap(data*factor,axis=0)/factor

# for traj in data.T[1:2]:
#     plt.plot(traj)

for traj in datauw.T:
    plt.plot(traj)

#plt.plot(np.average(data,axis=1))    
plt.show()
