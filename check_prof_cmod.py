#! /usr/bin/python

from pylab import *
from sys import argv,exit,stdout
import optparse as op
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline as US
from scipy import interpolate
import numpy as np
import re
from interp import *
from finite_differences import *
from read_EFIT_file_cmod import *
from read_cmod_pfile import *
from calc_gammaE import *
from read_cxrs_file import *
from read_iterdb_file import *
from w_iterdb import *
from matplotlib.backends.backend_pdf import PdfPages

e = 1.6*10**(-19)
a = 2.6863997038399379E-01
Zeff = 2.8
Zave = 10.

# flag: -o (output iterdb file)
parser = op.OptionParser()
parser.add_option('--out_iterdb','-o',action='store_const',const=1,default=0)
options,args = parser.parse_args()
out_iterdb =options.out_iterdb

efit_file_name = 'g1120907032.01012'
p_file_name = 'p1120907032.01012'
cxrs_file_name = 'cxrs1120907032.01010_v20140623'

# read EFIT file
psip_n, Rgrid, Zgrid, F, p, ffprime, pprime, psirz, qpsi, rmag, zmag, nw, psiax, psisep = read_EFIT_file_cmod(efit_file_name)
# calculate psip_n vs rhot_n
rho_tor_spl, rhot_n = calc_rho_tor(psip_n, psiax, psisep, qpsi, nw)
# calculate B_pol and B_tor at outboard midplane
psip_n_obmp, R_obmp, B_pol, B_tor = calc_B_fields(Rgrid, rmag, Zgrid, zmag, psirz, psiax, psisep, F,nw,psip_n)

psip_obmp = psip_n_obmp*(psisep-psiax)+psiax
rhot_n_obmp = rho_tor_spl(psip_obmp)

### read from p-file
### ne, ni are in the unit of 10^20 m^(-3)
### te, ti are in the unit of KeV
psi0, ne, te, ni, ti_pfile, nz_pfile=read_cmod_pfile(p_file_name,shift_Ti=False,shift_Te=False,add_impurity=True)

# convert into SI units
ne = ne*1E20
te = te*e*1E03
#ni = ni*1E20
#ti_pfile = ti_pfile*e*1E03
#nz_pfile = nz_pfile*1E20
# find ni,ti,ne,te,q on the grid of uniform psip_n at outboard midplane
rhot0=rho_tor_spl(psi0*(psisep-psiax)+psiax)

#ni_obmp = interp(psi0,ni,psip_n_obmp)
#ti_obmp = interp(psi0,ti,psip_n_obmp)
ne_obmp = interp(psi0,ne,psip_n_obmp)
te_obmp = interp(psi0,te,psip_n_obmp)

rhotp_obmp = interp(psi0,rhot0,psip_n_obmp)
q_obmp = interp(psip_n, qpsi, psip_n_obmp)

# shift cxrs measurements out by 0.005 psip_n
cxrs_shift = 0.005  #Shift in poloidal flux coord
psip_tb,tb,psip_nb,nb,psip_Er,Er = read_cxrs_file(cxrs_file_name)
psip_tb = psip_tb+cxrs_shift
#psip_nb = psip_nb+cxrs_shift
psip_Er = psip_Er+cxrs_shift
# convert into SI units
tb = tb*e*1E03
#nb = nb*1E18
Er = Er*1E03

tb_full = interp(psip_tb,tb,psi0)
Erplots=PdfPages('cmod_profiles.pdf')
#plt.plot(rhot0,tb_full/e*1E-3,'.',label='Ti(T_Boron shifted)')
#plt.plot(rhot0,te/e*1E-3,'.',label='Te')
plt.plot(psi0,tb_full/e*1E-3,'+-',label='Ti (T_Boron cxrs_shifted)')
plt.plot(psi0,te/e*1E-3,'+-',label='Te (p-file)')
plt.xlabel('psi_pol_n')
plt.ylabel('keV')
plt.legend()
plt.grid()
plt.axis((0.95,1.,0.,1.))
#plt.show()
Erplots.savefig()
plt.close()


ni_new = ne*(Zave-Zeff)/(Zave-1.)
nz_new = ne*(Zeff-1.)/Zave/(Zave-1.)
#plt.plot(rhot0,ne*1E-20,'.',label='ne')
#plt.plot(rhot0,ni_new*1E-20,'.',label='ni')
#plt.plot(rhot0,nz_new*1E-20,'.',label='nz')
plt.plot(psi0,ne*1E-20,'+-',label='ne (p-file)')
plt.plot(psi0,ni_new*1E-20,'+-',label='ni (ne*0.8)')
plt.plot(psi0,nz_new*1E-20,'+-',label='nz (ne*0.02)')
plt.xlabel('psi_pol_n')
plt.ylabel('10^20 m^-3')
plt.legend()
plt.grid()
plt.axis((0.95,1.,0.,1.5))
#plt.show()
Erplots.savefig()
plt.close()

psi_Er_sep = np.argmin(abs(psip_Er-1.))
psi_ped = psip_Er[:psi_Er_sep+1].copy()
Er_ped = Er[:psi_Er_sep+1].copy()
plt.plot(psi_ped,Er_ped,'+-',label='Er (cxrs_shifted)')
plt.xlabel('psi_pol_n')
plt.ylabel('kV/m')
plt.legend(loc=2)
#plt.show()
Erplots.savefig()
plt.close()
Erplots.close()

psi_Er_sep = np.argmin(abs(psip_Er-1.))
psi_ped = psip_Er[:psi_Er_sep].copy()
Er_ped = Er[:psi_Er_sep].copy()
R_ped = interp(psip_n_obmp,R_obmp,psi_ped)
Bp_ped = interp(psip_n_obmp,B_pol,psi_ped)
omega_tor = Er_ped/R_ped/Bp_ped
rhot_ped=rho_tor_spl(psi_ped*(psisep-psiax)+psiax)
#plt.plot(rhot_ped,omega_tor,label='new')
#plt.legend()
#plt.show()


omega_tor_full = interp(rhot_ped,omega_tor,rhot0)
