#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from sys import argv,exit,stdout
import optparse as op
from ParIO import *
import os
import sys
import optparse as op
from subprocess import call
from read_write_geometry import *
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

parser = op.OptionParser()
options,args = parser.parse_args()
geom_file = args[0]
gpars,geometry = read_geometry_global(geom_file)

nz,nx=np.shape(geometry['gxx'])
#geom_list = ["gxx","gxy","gxz","gyy","gyz","gzz","Bfield","dBdx","dBdy","dBdz","jacobian","geo_R","geo_Z","geo_c1","geo_c2"]
geom_list = ['dBdx','dBdz']
for i in geom_list:
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    z=np.arange(nz)
    x=np.arange(nx)
    z,x = np.meshgrid(z,x)
    surf = ax.plot_surface(z,x,geometry[i],rstride=1,cstride=1,cmap=cm.coolwarm,linewidth=0,antialiased=False)
    fig.colorbar(surf,shrink=0.5,aspect=5)
    plt.title(i)
    plt.show()

#avg_geom = {}
#for i in geom_list :
#    geom = geometry[i]
#    avg_geom[i] = np.empty(0)
#    for j in np.arange(nx):
#        avg = np.mean(geom[:,j])
#        avg_geom[i] = np.append(avg_geom[i],avg)
#    plt.plot(avg_geom[i])
#    plt.title(i)
#    plt.show()
