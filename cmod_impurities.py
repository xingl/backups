#! /usr/bin/python

from pylab import *
from sys import argv,exit,stdout
import optparse as op
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline as US
from scipy import interpolate
import numpy as np
import re
#from interp import *
#from finite_differences import *
#from read_EFIT_file_cmod import *
from read_cmod_pfile import *
#from calc_gammaE import *
from read_cxrs_file import *
#from read_iterdb_file import *
#from w_iterdb import *
#from matplotlib.backends.backend_pdf import PdfPages

e = 1.6*10**(-19)
a = 2.6863997038399379E-01
Zeff=2.8
Zave=10.

# flag: -o (output iterdb file)
#parser = op.OptionParser()
#parser.add_option('--out_iterdb','-o',action='store_const',const=1,default=0)
#options,args = parser.parse_args()
#out_iterdb =options.out_iterdb

efit_file_name = 'g1120907032.01012'
p_file_name = 'p1120907032.01012'
cxrs_file_name = 'cxrs1120907032.01010_v20140623'

### read from p-file
### ne, ni are in the unit of 10^20 m^(-3)
### te, ti are in the unit of KeV
psi0, ne, te, ni, ti_pfile, nz_pfile=read_cmod_pfile(p_file_name,shift_Ti=False,shift_Te=False,add_impurity=True)
#psipne, ne, psipte, te, psipni, ni, psipti, ti=read_cmod_pfile_raw(p_file_name)

# shift cxrs measurements out by 0.005 psip_n
cxrs_shift = 0.005  #Shift in poloidal flux coord
psip_tb,tb,psip_nb,nb,psip_Er,Er = read_cxrs_file(cxrs_file_name)
psip_tb = psip_tb+cxrs_shift
psip_nb = psip_nb+cxrs_shift
psip_Er = psip_Er+cxrs_shift

##plt.plot(psi0,ne,label='ne_pfile')
##plt.plot(psi0,ni,label='ni_pfile')
##plt.plot(psi0,nz_pfile,label='nz_pfile')
###plt.plot(psip_nb,nb*1E-02,label='shifted nb_cxrs')
###plt.plot(psip_nb,2*nb*1E-02,label='n_mo')
##plt.plot(psip_nb,3.*nb*1E-02,label='n_mo+n_b')
##plt.ylabel('10^20 m^-3')
##plt.legend()
##plt.show()

sep_ind = np.argmin(abs(psip_nb-1.))
psi_ped = psip_nb[:sep_ind+1].copy()

ni_ped = interp(psi0,ni,psi_ped)
#plt.plot(psi0,ni,label='ni_pfile')
#plt.plot(psi_ped,ni_ped,label='ni_ped')
#plt.legend()
#plt.show()
ne_ped = interp(psi0,ne,psi_ped)
#plt.plot(psi0,ne,label='ne_pfile')
#plt.plot(psi_ped,ne_ped,label='ne_ped')
#plt.legend()
#plt.show()
nb_ped = interp(psip_nb,nb,psi_ped)
#print np.mean(nb_ped*1E-02/ne_ped)
#plt.plot(psip_nb,nb,label='ne_pfile')
#plt.plot(psi_ped,nb_ped,label='ne_ped')
#plt.legend()
#plt.show()

#nz1=(ne_ped-ni_ped)/Zave
#plt.plot(psi_ped,ni_ped/ne_ped,label='ni/ne pfile')
#plt.plot(psi_ped,nz1/ne_ped,label='nz/ne pfile')
#plt.show()
#print np.mean(ni_ped/ne_ped)
#print np.mean(nz1/ne_ped)

print 'ni/ne = ', (Zave-Zeff)/(Zave-1.)
print 'nz/ne = ', (Zeff-1.)/Zave/(Zave-1.)
ni_new = ne*(Zave-Zeff)/(Zave-1.)
plt.plot(psi0,ni,label='ni_pfile')
plt.plot(psi0,ni_new,label='ni_new')
plt.plot(psi0,ne,label='ne_pfile')
plt.legend()
plt.show()
