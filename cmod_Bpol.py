#! /usr/bin/python

from pylab import *
from sys import argv,exit,stdout
import optparse as op
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline as US
from scipy import interpolate
import numpy as np
#import re
#from interp import *
#from finite_differences_x import *
from read_EFIT_file_cmod import *
#from read_EFIT_file_cmod_backup import *
#from read_cmod_pfile import *
#from calc_gammaE import *
#from read_cxrs_file import *
#from read_iterdb_file import *
#from w_iterdb import *
#from matplotlib.backends.backend_pdf import PdfPages

e = 1.6*10**(-19)
a = 2.6863997038399379E-01

# flag: -o (output iterdb file)
parser = op.OptionParser()
#parser.add_option('--out_iterdb','-o',action='store_const',const=1,default=0)
options,args = parser.parse_args()
#out_iterdb =options.out_iterdb

efit_file_name = 'g1120907032.01012'
p_file_name = 'p1120907032.01012'
cxrs_file_name = 'cxrs1120907032.01010_v20140623'

# read EFIT file
psip_n, Rgrid, Zgrid, F, p, ffprime, pprime, psirz, qpsi, rmag, zmag, nw, psiax, psisep = read_EFIT_file_cmod(efit_file_name)
plt.plot(psip_n,qpsi,label='q')
plt.legend()
plt.axis((0.95,1.,3.,5.))
plt.show()

# calculate psip_n vs rhot_n
rho_tor_spl, rhot_n = calc_rho_tor(psip_n, psiax, psisep, qpsi, nw)
# calculate B_pol and B_tor at outboard midplane
psip_n_obmp, R_obmp, B_pol, B_tor = calc_B_fields(Rgrid, rmag, Zgrid, zmag, psirz, psiax, psisep, F,nw,psip_n)
#psip_n_obmp2, R_obmp2, B_pol2, B_tor2 = calc_B_fields_backup(Rgrid, rmag, Zgrid, zmag, psirz, psiax, psisep, F,nw,psip_n)

psip_obmp = psip_n_obmp*(psisep-psiax)+psiax
rhot_n_obmp = rho_tor_spl(psip_obmp)

plt.plot(psip_n_obmp,B_pol,label='Bpol')
#plt.plot(psip_n_obmp2,B_pol2,label='Bpol2')
plt.axis((0.95,1.0,0.9,1.1))
plt.legend()
plt.show()
plt.plot(psip_n_obmp,B_tor,label='Btor')
#plt.plot(psip_n_obmp2,B_tor2,label='Btor2')
plt.axis((0.95,1.,4.,5.))
plt.legend()
plt.show()

ped_t_ind = np.argmin(abs(rhot_n-0.95))
sep_ind = np.argmin(abs(rhot_n-0.99))
rhot_n_ped = rhot_n[ped_t_ind:sep_ind+1].copy()
q_ped = qpsi[ped_t_ind:sep_ind+1].copy()
q_spl = interpolate.UnivariateSpline(rhot_n_ped,q_ped)
shat = q_spl(rhot_n_ped,nu=1)*rhot_n_ped/q_spl(rhot_n_ped)
plt.plot(rhot_n,qpsi,label='before')
plt.plot(rhot_n_ped,q_spl(rhot_n_ped),label='after')
plt.axis((0.95,1.,4.,5.))
plt.legend()
plt.show()
plt.plot(rhot_n_ped,shat,label='shat')
plt.xlabel('rhot_n')
plt.legend()
plt.show()
