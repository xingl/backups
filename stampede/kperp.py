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
ly = float(args[1])
kymin = 2*np.pi/ly
lx = float(args[2])
dkx = 2*np.pi/lx
gpars,geometry = read_geometry_local(geom_file)

nz=np.shape(geometry['ggxx'])
z_grid_center = np.linspace(-1.,1.,nz[0],endpoint=False)
ggxx = geometry['ggxx']
ggxy = geometry['ggxy']
ggyy = geometry['ggyy']

#for i in np.arange(-5,6):
for i in np.arange(-10,11):
	kx = i*dkx
	print 'kx = ', kx
	print 'kymin = ', kymin
	#k1 = (kx*ggxx+kymin*ggxy)/np.sqrt(ggxx)
	#k2 = np.sqrt((ggxx*ggyy-ggxy**2)/ggxx)*kymin
	#kperp = np.sqrt(k1**2+k2**2)
	kperp2 = np.sqrt(ggxx*kx**2+2.*ggxy*kx*kymin+ggyy*kymin**2)
	#plt.plot(z_grid,kperp)
	z_grid = z_grid_center + i*2.
	plt.plot(z_grid,kperp2)
plt.axvline(x=0.,color='r')
plt.xlabel('z/pi')
plt.ylabel('|k_perp|')
plt.title('kymin = '+str(kymin))
plt.grid()
plt.show()


