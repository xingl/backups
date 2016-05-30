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
from finite_differences_x import *
from read_EFIT_file_cmod import *
from read_cmod_pfile import *
from calc_gammaE import *
from read_cxrs_file import *
from read_iterdb_file import *
#from w_iterdb import *
from write_iterdb_new import *
from matplotlib.backends.backend_pdf import PdfPages

M = 3.3*10**(-27)
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
###ne = ne*1E20
###te = te*1E03
#ni = ni*1E20
#ti_pfile = ti_pfile*e*1E03
#nz_pfile = nz_pfile*1E20
# find ni,ti,ne,te,q on the grid of uniform psip_n at outboard midplane
rhot0=rho_tor_spl(psi0*(psisep-psiax)+psiax)

#ni_obmp = interp(psi0,ni,psip_n_obmp)
#ti_obmp = interp(psi0,ti,psip_n_obmp)
ne_obmp = interp(psi0,ne,psip_n_obmp)
te_obmp = interp(psi0,te,psip_n_obmp)

lamD_rhoe = sqrt((B_pol**2+B_tor**2)/(ne_obmp*10))
plt.plot(psip_n_obmp,lamD_rhoe,label='lambda_D/rho_e')
plt.legend(loc=2)
plt.axis((0.95,1.,1.,2.))
plt.show()
#rhotp_obmp = interp(psi0,rhot0,psip_n_obmp)
#q_obmp = interp(psip_n, qpsi, psip_n_obmp)

# shift cxrs measurements out by 0.005 psip_n
cxrs_shift = 0.005  #Shift in poloidal flux coord
psip_tb,tb,psip_nb,nb,psip_Er,Er = read_cxrs_file(cxrs_file_name)
psip_tb = psip_tb+cxrs_shift
#psip_nb = psip_nb+cxrs_shift
psip_Er = psip_Er+cxrs_shift
Er = Er*1E03

tb_full = interp(psip_tb,tb,psi0)
plt.plot(rhot0,tb_full,'+-',label='Ti(T_Boron shifted)')
plt.plot(rhot0,te,'+-',label='Te')
plt.xlabel('rhot')
plt.ylabel('keV')
plt.legend()
plt.grid()
plt.axis((0.95,1.,0.,1.))
plt.show()

plt.plot(rhot0,te/tb_full,'+-',label='Te/Ti')
plt.xlabel('rhot')
plt.legend()
plt.grid()
plt.axis((0.95,1.,0.,1.5))
plt.show()

ni_new = ne*(Zave-Zeff)/(Zave-1.)
nz_new = ne*(Zeff-1.)/Zave/(Zave-1.)
plt.plot(rhot0,ne,'+-',label='ne')
plt.plot(rhot0,ni_new,'+-',label='ni')
plt.plot(rhot0,nz_new,'+-',label='nz')
plt.xlabel('rhot')
plt.ylabel('10^20 m^-3')
plt.legend()
plt.grid()
plt.axis((0.95,1.,0.,1.5))
plt.show()


psi_Er_sep = np.argmin(abs(psip_Er-1.))
psi_Er_top = np.argmin(abs(psip_Er-0.96))
psi_ped = psip_Er[psi_Er_top:psi_Er_sep+1].copy()
Er_ped = Er[psi_Er_top:psi_Er_sep+1].copy()
plt.plot(psi_ped,Er_ped,'+-',label='Er_cxrs')
plt.ylabel('V/m')
plt.legend()
plt.legend()
plt.show()

R_ped = interp(psip_n_obmp,R_obmp,psi_ped)
Bp_ped = interp(psip_n_obmp,B_pol,psi_ped)
omega_tor = Er_ped/R_ped/Bp_ped
rhot_ped=rho_tor_spl(psi_ped*(psisep-psiax)+psiax)

omega_tor_full = interp(rhot_ped,omega_tor,rhot0)

#plt.plot(rhot0,omega_tor_full,'+-',label='omega_tor')
#plt.plot(rhot_ped,omega_tor,'+-',label='omega_tor_ped')
#plt.legend(loc=2)
#plt.grid()
#plt.axis((0.95,1.,-1.5E5,100.))
#plt.show()

rhot1 = np.linspace(rhot_ped[0],rhot_ped[-1],300)
omega_tor_spl = US(rhot_ped,omega_tor)
q_ped = interp(psip_n,qpsi,psi_ped)
te_ped = interp(psi0,te,psi_ped)
q_rhot_spl = US(rhot_ped,q_ped)
te_rhot_spl = US(rhot_ped,te_ped)

omega_t = np.sqrt(te_rhot_spl(rhot1)*1E3*e/M)/a
plt.plot(rhot1,omega_t,'+-',label='omega_t')
plt.legend()
plt.grid()
plt.show()

Er_rhot_spl = US(rhot_ped,Er_ped)
Bt_ped = interp(psip_n_obmp,B_tor,psi_ped)
B_ped = np.sqrt(Bt_ped**2+Bp_ped**2)
B_rhot_spl = US(rhot_ped,B_ped)
plt.plot(rhot_ped,B_ped,label='before')
plt.plot(rhot1,B_rhot_spl(rhot1),label='after')
plt.legend()
plt.show()
kymin=0.1
omega_tor1 = kymin*Er_rhot_spl(rhot1)/B_rhot_spl(rhot1)/omega_t
plt.plot(rhot1,omega_tor1,'+-',label='omega_tor1')
plt.legend(loc=2)
plt.grid()
plt.show()

gamE = rhot1/q_rhot_spl(rhot1)*omega_tor_spl(rhot1,nu=1)/omega_t

plt.plot(rhot1,abs(gamE),'+-',label='ExB rate')
plt.legend(loc=2)
plt.grid()
plt.show()

if out_iterdb:
    #f = open('Er_profile.dat','w')
    #f.write('#1.rhot_n 2.gammaE\n')
    #np.savetxt(f,np.column_stack((rhot_gene,gammaE_gene)))

    #f = open('omega_profile.dat','w')
    #f.write('#1.rhot_n 2.psip_n, 3.gammaE\n')
    #np.savetxt(f,np.column_stack((rhot0,psi0,omega_tor_full)))
    #f.close()
    
    #file_out_base = 'profiles_wb' 
    base_number = '1120907032.01012'
#    output_iterdb(rhot0,ne,te,ni,ti/e,file_out_base+base_number,base_number,'9999',vrot=omega_tor_full,nimp=nz)
    #write_iterdb routine needs temperature in KeV and density in 10^19 m^-3
    output_iterdb(rhot0,psi0,ne*10,te,ni_new*10,tb_full,'cmod'+base_number,base_number,'9999',vrot=omega_tor_full,nimp=nz_new*10)

    #rhop=np.sqrt(psi0)
#    f = open('gene_profiles'+file_out_base+'_gene_e.dat','w')
    #f = open('profiles_e','w')
    #f.write('# 1.rho_tor 2.rho_pol 3.Te(kev) 4.ne(10^19m^-3)\n#\n')
    #np.savetxt(f,np.column_stack((rhot0,rhop,te*1E-03,ne*1E-19)))
    #f.close()

#    f = open('gene_profiles'+file_out_base+'_gene_i.dat','w')
    #f = open('profiles_i','w')
    #f.write('# 1.rho_tor 2.rho_pol 3.Ti(kev) 4.ni(10^19m^-3)\n#\n')
#    np.savetxt(f,np.column_stack((rhot0,psi0,ti/e,ni)))
    #np.savetxt(f,np.column_stack((rhot0,rhop,tb_full*1E-03,ni*1E-19)))
    #f.close()
		
    #f = open('profiles_z','w')
    #f.write('# 1.rho_tor 2.rho_pol 3.Ti(kev) 4.nz(10^19m^-3)\n#\n')
#    np.savetxt(f,np.column_stack((rhot0,psi0,ti/e,ni)))
    #np.savetxt(f,np.column_stack((rhot0,rhop,tb_full*1E-03,nz*1E-19)))
    #f.close()
