import numpy as np
import os
import sys

"""
Helper script for tiling the glass network to get larger system sizes
Usage: tilexy.py doublings_along_x doublings_along_y
"""

ZEROCORR=0

def bondlengths(coords,bonds):
    return np.array([np.linalg.norm(coords[i]-coords[j]) for i,j in bonds])

def cross_bdry(coords,bonds,boxsize,axis=0):
    ''' find bonds that cross axis boundary '''
    deltax = (coords[bonds[:,1]]-coords[bonds[:,0]])[:,axis]
    idxslr = np.where(deltax > boxsize/2.)
    idxsrl = np.where(deltax < -boxsize/2.)
    return idxslr, idxsrl

    
def find_edges(coords,bonds,dx,dy):
    idxs_x1, idxs_x2 = cross_bdry(coords,bonds,dx,axis=0)
    idxs_y1, idxs_y2 = cross_bdry(coords,bonds,dy,axis=1)
    
    return idxs_x1, idxs_x2, idxs_y1, idxs_y2
    

def draw_bonds(coords, bonds,scale_factor=.3,resolution=6, cutoff=2):
    '''
    scale_factor sets size of the spheres drawn at each point
    resolution sets the detail of the spheres drawn at each point. higher = slower code
    '''
    x,y = coords.T
    z = np.zeros_like(x)  # to plot 2d points, include z coordinate of 0
    points = mlab.points3d(x,y,z, scale_mode='none',scale_factor=scale_factor,resolution=resolution)
    
    
    # eliminate really long bonds because of periodic boundary conditions
    if cutoff is not None:
        bl = np.array(bondlengths(coords,bonds))
        bonds =  bonds[np.where(bl < cutoff)]
    
    points.mlab_source.dataset.lines = bonds
    points.mlab_source.update()
    mlab.pipeline.surface(points, representation='wireframe',color=(0.5,0.5,0.5))
    
def write_bonds(bds, bdsfull, bdname):
    np.savetxt('__bds', bds+ZEROCORR, fmt='%d\t%d')
    np.savetxt('__bddata', bdsfull[:,2:], fmt="%d\t%d\t%2.16f")
    os.system("paste __bds __bddata > %s" % bdname)
    os.system("rm -f __bds __bddata")
    

pts = np.loadtxt('glass.pts')
bdsfull = np.loadtxt('glass.bds')
bds = np.array(bdsfull[:,:2]-ZEROCORR,dtype=int)
dx = 76.68
dy = 59.0282

n0 = len(pts)
m0 = len(bds)
ntilex = 1
ntiley = 1
if (len(sys.argv) == 3):
    ntilex = int(sys.argv[1])
    ntiley = int(sys.argv[2])

ptname = 'glass_tile%d_%d.pts' % (ntilex,ntiley)
bdname = 'glass_tile%d_%d.bds' % (ntilex,ntiley)
wallbdname = 'glass_tile%d_%d_wall.bds' % (ntilex,ntiley) # if using wall move rather than box shear, this creates a mesh with free ends along the x direction
wallname = 'wall_tile%d_%d.txt' % (ntilex,ntiley) # if using wall move rather than box shear, this identifies the points on the left edge wwhich must be sheared to measure propagation

for i in range(ntilex):
    shift = np.zeros_like(pts)
    shift[:,0] = 1
    
    idxs_x1, idxs_x2, idxs_y1, idxs_y2 = find_edges(pts,bds,dx,dy)
    
    pts = np.concatenate((pts, pts+shift*dx),axis=0)
    bds = np.concatenate((bds, bds+n0), axis=0)
    bdsfull = np.concatenate((bdsfull,bdsfull), axis=0) # update the bonds part at end
        
    # stitch the middle
    bds[idxs_x1] += np.array([n0,0])
    bds[idxs_x2] += np.array([0,n0])
    
    # stitch the edges
    bds[idxs_x1[0]+m0] -= np.array([n0,0])
    bds[idxs_x2[0]+m0] -= np.array([0,n0])
    
    dx *= 2
    n0 *= 2
    m0 *= 2
    
for i in range(ntiley):
    shift = np.zeros_like(pts)
    shift[:,1] = 1
    
    idxs_x1, idxs_x2, idxs_y1, idxs_y2 = find_edges(pts,bds,dx,dy)
    
    pts = np.concatenate((pts, pts+shift*dy),axis=0)
    bds = np.concatenate((bds, bds+n0), axis=0)
    bdsfull = np.concatenate((bdsfull,bdsfull), axis=0) # update the bonds part at end
        
    # stitch the middle
    bds[idxs_y1] += np.array([n0,0])
    bds[idxs_y2] += np.array([0,n0])
    
    # stitch the edges
    bds[idxs_y1[0]+m0] -= np.array([n0,0])
    bds[idxs_y2[0]+m0] -= np.array([0,n0])
    
    dy *= 2
    n0 *= 2
    m0 *= 2

# write bonds
write_bonds(bds,bdsfull,bdname)

np.savetxt(ptname,pts,fmt="%2.16f")

# create wall info on left edge for shear
idxs_x1, idxs_x2, idxs_y1, idxs_y2 = find_edges(pts,bds,dx,dy)
masses = np.ones(len(pts)) 
walls = np.zeros(len(pts),dtype=int) 
walls[bds[idxs_x1][:,0]] = 1
walls[bds[idxs_x2][:,1]] = 1

masswall = np.vstack((walls, masses))
np.savetxt(wallname,masswall.T, fmt ="%d\t%f") 

toclip = np.concatenate((idxs_x1[0],idxs_x2[0]))
print len(toclip)
bdsclip = np.delete(bds, toclip, axis=0)
bdsfullclip = np.delete(bdsfull, toclip, axis=0)
write_bonds(bdsclip,bdsfullclip,wallbdname)

