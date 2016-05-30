from pylab import *
from sys import argv,exit,stdout
import optparse as op
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline as US
from scipy import interpolate
import numpy as np
from finite_differences import *
from interp import *
from read_EFIT_file_new import *
from fields_along_fluxsurface import *
#from read_rbsProfs import *
from calc_fields_from_EFIT_2 import *

parser = op.OptionParser()
options,args = parser.parse_args()
efit_file_name = args[0]
#rbsProfs_file_name = args[1]
rhot_fs = args[1]
psip_fs = rhot_fs
#psip_fs = rhot2psip(rbsProfs_file_name,rhot_fs)
#print psip_fs

R_fs, Z_fs, Bp_R, Bp_Z, B_pol, B_tor, B_tot = flux_surface_B_fields(efit_file_name,psip_fs)

plt.plot(R_fs,Z_fs,label='flux surface')
plt.xlabel('R')
plt.ylabel('Z')
plt.legend(loc=2)
plt.show()
plt.plot(R_fs,Bp_R,label='Bp_R')
plt.plot(R_fs,Bp_Z,label='Bp_Z')
plt.xlabel('R')
plt.ylabel('B')
plt.legend(loc=2)
plt.show()

curr = 0.
for i in np.arange(len(R_fs)-1):
    this_curr = 0.5*(Bp_R[i]+Bp_R[i+1])*(R_fs[i+1]-R_fs[i])\
		+0.5*(Bp_Z[i]+Bp_Z[i+1])*(Z_fs[i+1]-Z_fs[i])
    curr = curr + this_curr/4./np.pi*1E7
print 'I(MA) = ', curr*1E-6
