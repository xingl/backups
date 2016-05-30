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
kymin = float(args[1])

gpars,geometry = read_geometry_local(geom_file)


q =  gpars['q0']
Lref = gpars['Lref']
cy = gpars['Cy']
gyy = geometry['ggyy']
gyz = geometry['ggyz']
gzz = geometry['ggzz']

#kymin =   0.5000E-01
rhostar  =   0.16839125E-02
print Lref*rhostar
k_theta = q*cy/np.sqrt(q**2*cy**2*gyy+2*q*cy*gyz+gzz)*kymin/Lref/rhostar
zgrid = np.linspace(-1,1,int(gpars['gridpoints']),endpoint=False)
#plt.plot(zgrid,k_theta)
#plt.show()
print 'k_theta(cm^-1) =', k_theta[int(gpars['gridpoints'])/2]/100
