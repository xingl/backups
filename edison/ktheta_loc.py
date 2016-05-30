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

M = 2.*1.6726219E-27
c = 3.E08
e = 1.6021765E-19
parser = op.OptionParser()
options,args = parser.parse_args()
suffix = '_'+args[0]
geom_file = 'tracer_efit'+suffix



gpars,geometry = read_geometry_local(geom_file)

par = Parameters()
par.Read_Pars('parameters'+suffix)
pars = par.pardict

kymin = pars['kymin']

q =  gpars['q0']
Lref = gpars['Lref']
cy = gpars['Cy']
gyy = geometry['ggyy']
gyz = geometry['ggyz']
gzz = geometry['ggzz']

#kymin =   0.5000E-01
#rhostar  =   0.13735059E-02#0.16839125E-02
rhostar = pars['rhostar']
print 'Lref*rhostar = ',Lref*rhostar
print 'sqrt(Te/M)*M/e/B = ',np.sqrt(pars['Tref']*1E03*e/M)*M/e/pars['Bref']
k_theta = q*cy/np.sqrt(q**2*cy**2*gyy+2*q*cy*gyz+gzz)*kymin/Lref/rhostar
zgrid = np.linspace(-1,1,int(gpars['gridpoints']),endpoint=False)
#plt.plot(zgrid,k_theta)
#plt.show()
print 'k_theta(cm^-1) =', k_theta[int(gpars['gridpoints'])/2]/100
print kymin
n0 = int(kymin*cy/rhostar)
print 'int(kymin*cy/rhostar) = ',n0
print 'n0_global', pars['n0_global']
