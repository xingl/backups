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
ky = 2*np.pi/ly
lx = float(args[2])
dkx = 2*np.pi/lx

gpars,geometry = read_geometry_local(geom_file)

nz=np.shape(geometry['ggxx'])
z_grid_center = np.linspace(-1.,1.,nz[0],endpoint=False)
ggxy = geometry['ggxy']
ggxx = geometry['ggxx']
ggxy = geometry['ggxy']
ggxz = geometry['ggxz']
ggyy = geometry['ggyy']
ggyz = geometry['ggyz']
gdBdx = geometry['gdBdx']
gdBdy = geometry['gdBdy']
gdBdz = geometry['gdBdz']
gBfield = geometry['gBfield']

#for i in np.arange(-5,6):
for i in np.arange(-10,11):
	kx = dkx*i
	print 'kx = ', kx
	print 'ky = ', ky
	gamma1 = ggxx*ggyy-ggxy**2
	gamma2 = ggxx*ggyz-ggxy*ggxz
	gamma3 = ggxy*ggyz-ggyy*ggxz
	Kx = -(gdBdy+gamma2/gamma1*gdBdz)
	Ky = gdBdx-gamma3/gamma1*gdBdz

	#plt.plot(z_grid,Ky)
	#plt.ylabel('Ky')
	#plt.show()

	omega_d_y = Ky*ky
	omega_d_x = Kx*kx

	omega_d = omega_d_y + omega_d_x
	z_grid = z_grid_center + i*2

	plt.plot(z_grid,omega_d)
plt.axvline(x=0.,color='r')
plt.xlabel('z/pi')
plt.ylabel('omega_d')
plt.title('kymin = '+str(ky))
plt.grid()
plt.show()

