#!/usr/bin/env python



from __future__ import division

import numpy as np
from scipy.special import gamma
import itertools as itt
import sys

def get_boxes(nbins_lin):
    filler = np.frompyfunc(lambda x: list(), 1, 1)
    boxes = np.empty((nbins_lin*nbins_lin), dtype=np.object)
    filler(boxes, boxes)
    
    return boxes

def binit(pos, L, radius):
    nbins_lin = int(L/(2*radius))
    xbins = np.linspace(0,L,num=nbins_lin+1)
    ybins = xbins

    digposx = np.digitize(pos[:,0],xbins)-1
    digposy = np.digitize(pos[:,1],ybins)-1
    digpos = nbins_lin*digposy+digposx
    boxes=get_boxes(nbins_lin)

    for i in range(len(pos)):
        boxes[digpos[i]].append(i)

    return boxes, digpos, nbins_lin

def get_neighbor_boxes(digpos, nbins_lin):
# http://stackoverflow.com/questions/15679719/python-neighbors-on-a-regular-grid
    idy, idx = np.unravel_index(digpos, dims=(nbins_lin, nbins_lin))
    
    neigh_idx = np.vstack((idx-1, idx, idx+1, idx-1, idx, idx+1, idx-1, idx, idx+1)) # !!! 2d
    neigh_idy = np.vstack((idy-1, idy-1, idy-1, idy, idy, idy, idy+1, idy+1, idy+1))
    bindex = np.ravel_multi_index((neigh_idy, neigh_idx), dims=(nbins_lin,nbins_lin), mode='wrap').T
    return bindex


def get_particles_from_box(box_list):
    return [list(itt.chain(*[y for y in x])) for x in box_list]

def make_pairs(npart):
    j = np.array(list(itt.chain(*npart)))
    i = np.array(list(itt.chain(*[[i]*len(npart[i]) for i in range(len(npart))])))
    pairs = np.vstack((i,j)).T
    pairs = pairs[pairs[:,1]-pairs[:,0] != 0] # remove self pairs

    return pairs
# def get_neighbor_particles(digpos, nbins_lin):
    

def apply_PBC(separation, L, d):
    Lhalf = L/2
    
    for s in range(d):
        separation[separation[:,s]>Lhalf,s] -= L
        separation[separation[:,s]<-Lhalf,s] += L
    
def get_particles_to_move(pos, L, radius, d):
    boxes, digpos, nbins_lin = binit(pos, L, radius)

    nbox = get_neighbor_boxes(digpos,nbins_lin)
    npart = get_particles_from_box(boxes[nbox])

    pairs = make_pairs(npart)

    separations = pos[pairs[:,1]]-pos[pairs[:,0]]
    apply_PBC(separations,L,d)

    distances = np.linalg.norm(separations, axis=1)
    to_be_moved = pairs[distances<2*radius][:,0]
    return to_be_moved




def rad(N,L,phi,d):
    V = np.power(L,d)
    return np.power((gamma(d/2+1)/np.power(np.pi,d/2))*(V*phi/N), 1/d)

n = int(sys.argv[1])
dim = 2
l = 1
Phi = float(sys.argv[2])
radius = rad(n,l,Phi,dim)
positions = np.random.rand(n, dim)
to_be_moved = get_particles_to_move(positions, l, radius, dim)

iterations = 0
name = "phi"+str(Phi)+".dat"
out_file_overlaps = open("overlaps_"+name, "w")

while len(to_be_moved)>0:
    print(iterations,len(to_be_moved), radius)
    to_be_moved = get_particles_to_move(positions, l, radius, dim)
    positions[to_be_moved] = np.random.rand(len(to_be_moved), dim)
    iterations += 1

    out_file_overlaps.write(str(len(to_be_moved)/n)+"\n")
    if iterations%100==0:
        np.savetxt("positions_"+name, positions)

out_file_overlaps.close()
np.savetxt("positions_"+name, positions)

