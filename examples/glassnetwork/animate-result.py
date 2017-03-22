#!/usr/bin/env python2

import numpy as np
import pylab
import sys
import os
import glob

pylab.ion()

def bondlengths(coords,bonds):
    return np.array([np.linalg.norm(coords[i]-coords[j]) for i,j in bonds])

    
def draw_mesh(filename):
    if filename[-1] == '.': filename += 'pts'
    print filename
    ar = np.loadtxt(filename)
    print ar.shape
    
    bondfile = filename[:-4] + '.bds'
    connections = np.loadtxt(bondfile)
    bonds = np.array(connections[:,:2],dtype=int)  # M x 2 array, discard remaining columns in the connection file
    bonds = bonds-np.min(bonds)  # index of points starts from 1 in the data files
    
    vis2d(ar,bonds)

def vis2d(pts,bds,draw_points=False,draw_edges=True,eigenfunction=None,eigenfunctions=None,ecolors=None,figname=None,noshow=False,figidx=None):
    fig = pylab.figure(figsize=(8,8))
    ax = fig.add_subplot(1, 1, 1, xticks=[], yticks=[])
    fig.figurePatch.set_edgecolor('white')
    fig.figurePatch.set_facecolor('white')
    ax.set_aspect('equal')
    ax.set_axis_off()
    
    
    if (draw_edges):
        cutoff = 10.
        xpairs = [[pts[bd[0]][0], pts[bd[1]][0]] for bd in bds if np.linalg.norm(pts[bd[0]]-pts[bd[1]])<cutoff]
        ypairs = [[pts[bd[0]][1], pts[bd[1]][1]] for bd in bds if np.linalg.norm(pts[bd[0]]-pts[bd[1]])<cutoff]
        xlist = []
        ylist = []
        for xends, yends in zip(xpairs,ypairs):
            xlist.extend(xends)
            xlist.append(None)
            ylist.extend(yends)
            ylist.append(None)
        plotpts, = pylab.plot(xlist,ylist,'-',color=(.5,.5,.5))
        
    if (draw_points):
        pylab.scatter(pts[:,0],pts[:,1])
        
    if (eigenfunction is not None):
        pylab.quiver(pts[:,0], pts[:,1], eigenfunction[::2], eigenfunction[1::2],color='r')
    
    if (eigenfunctions is not None):
        for i,eig in enumerate(eigenfunctions):
            pylab.quiver(pts[:,0], pts[:,1], eig[::2], eig[1::2],color=ecolors[i])
    
    if figname is None and not noshow:
        pylab.show()
        
    if figname is not None:
        pylab.savefig(figname)
        pylab.close()
        
    return plotpts

#############################################


coordfiletemplate = 'sim%s.pts'
bondfiletemplate = 'sim%s.bds'

frame = 0
if len(sys.argv) > 1: frame = int(sys.argv[1])

coordfile = coordfiletemplate % (frame)
bondfile = bondfiletemplate % (frame)

print coordfile,bondfile
# read point coordinates and bond connections
coords = np.loadtxt(coordfile) # N x 2 array
connections = np.loadtxt(bondfile)
bonds = np.array(connections[:,:2],dtype=int)  # M x 2 array, discard remaining columns in the connection file
bonds = bonds-np.min(bonds)  # index of points starts from 1 in the data files
cutoff = 10.
bl = np.array(bondlengths(coords,bonds))
bonds =  bonds[np.where(bl < cutoff)]
plotpts = vis2d(coords,bonds,noshow=True)
pylab.draw()

for i in range(len(glob.glob(coordfiletemplate % '*'))):
    coordfile = coordfiletemplate % i
    print coordfile
    pts = np.loadtxt(coordfile)
    bds = bonds
    
    xpairs = [[pts[bd[0]][0], pts[bd[1]][0]] for bd in bds if np.linalg.norm(pts[bd[0]]-pts[bd[1]])<cutoff]
    ypairs = [[pts[bd[0]][1], pts[bd[1]][1]] for bd in bds if np.linalg.norm(pts[bd[0]]-pts[bd[1]])<cutoff]
    xlist = []
    ylist = []
    for xends, yends in zip(xpairs,ypairs):
        xlist.extend(xends)
        xlist.append(None)
        ylist.extend(yends)
        ylist.append(None)
    
    plotpts.set_xdata(xlist)
    plotpts.set_ydata(ylist)
    pylab.draw()
    
#~ raw_input('press Enter')
