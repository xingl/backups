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
ni = ni*1E20
ti_pfile = ti_pfile*e*1E03
nz_pfile = nz_pfile*1E20
# find ni,ti,ne,te,q on the grid of uniform psip_n at outboard midplane
rhot0=rho_tor_spl(psi0*(psisep-psiax)+psiax)

ni_obmp = interp(psi0,ni,psip_n_obmp)
#ti_obmp = interp(psi0,ti,psip_n_obmp)
ne_obmp = interp(psi0,ne,psip_n_obmp)
te_obmp = interp(psi0,te,psip_n_obmp)

rhotp_obmp = interp(psi0,rhot0,psip_n_obmp)
q_obmp = interp(psip_n, qpsi, psip_n_obmp)

# shift cxrs measurements out by 0.005 psip_n
cxrs_shift = 0.005  #Shift in poloidal flux coord
psip_tb,tb,psip_nb,nb,psip_Er,Er = read_cxrs_file(cxrs_file_name)
psip_tb = psip_tb+cxrs_shift
psip_nb = psip_nb+cxrs_shift
psip_Er = psip_Er+cxrs_shift
# convert into SI units
tb = tb*e*1E03
nb = nb*1E18
Er = Er*1E03

tb_full = interp(psip_tb,tb,psi0)
plt.plot(rhot0,tb_full/e*1E-3,'.',label='Ti_shifted')
plt.plot(rhot0,te/e*1E-3,'.',label='Te')
plt.xlabel('rhot')
plt.ylabel('keV')
plt.legend()
plt.grid()
plt.axis((0.95,1.,0.,1.))
plt.show()

Zeff=2.8
Z=(ne*Zeff-ni)/(ne-ni)
nz=(ne-ni)/Z
plt.plot(rhot0,Z,label='Z_fromZeff=2.8')
(x1,x2,y1,y2)=plt.axis()
plt.axis((0.95,1.,y1,y2))
plt.legend()
plt.xlabel('rhot')
plt.show()
plt.plot(rhot0,nz*1E-20,'.',label='nz')
plt.xlabel('rhot')
plt.ylabel('10^20 m^-3')
plt.legend()
plt.grid()
plt.axis((0.95,1.,0.,0.02))
plt.show()
#nb_full = interp(psip_nb,nb,psi0)
plt.plot(rhot0,ni*1E-20,'.',label='ni')
plt.plot(rhot0,ne*1E-20,'.',label='ne')
plt.plot(rhot0,(ni+Z*nz)*1E-20,'.',label='ni+nz*Z')
plt.xlabel('rhot')
plt.ylabel('10^20 m^-3')
plt.legend()
plt.grid()
plt.axis((0.95,1.,0.,1.5))
plt.show()

psi_Er_f = np.argmin(abs(psip_n_obmp-psip_Er[0]))
if psip_Er[-1] <= psip_n_obmp[-1] :
    psi_Er_l = np.argmin(abs(psip_n_obmp-psip_Er[-1]))
else:
    psi_Er_l = np.argmin(abs(psip_n_obmp-1.01))

psipn_Er = np.linspace(psip_n_obmp[psi_Er_f],psip_n_obmp[psi_Er_l],500)
Er_obmp = interp(psip_Er,Er,psipn_Er)
rhotn_Er = rho_tor_spl(psipn_Er*(psisep-psiax)+psiax)

rhot_gene,gammaE_gene,omega_tor_gene = calc_gammaE_gene(R_obmp[psi_Er_f:psi_Er_l+1], rhot_n_obmp[psi_Er_f:psi_Er_l+1], te_obmp[psi_Er_f:psi_Er_l+1], B_pol[psi_Er_f:psi_Er_l+1], Er_obmp, rhotn_Er,q_obmp[psi_Er_f:psi_Er_l+1], a)
# convert omega that's valid on [0.905,1.01] to [0,1] so to write out
omega_tor_full = interp(rhot_gene,omega_tor_gene,rhot0)
#Erprof = np.genfromtxt('Er_profile.dat')

#plt.plot(rhot_gene,gammaE_gene,label='new')
#plt.plot(Erprof[:,0],Erprof[:,1],label='old')
#plt.show()


if out_iterdb:
    f = open('Er_profile.dat','w')
    f.write('#1.rhot_n 2.gammaE\n')
    np.savetxt(f,np.column_stack((rhot_gene,gammaE_gene)))

    f = open('omega_profile.dat','w')
    f.write('#1.rhot_n 2.psip_n, 3.gammaE\n')
    np.savetxt(f,np.column_stack((rhot0,psi0,omega_tor_full)))
    f.close()
    
    file_out_base = 'profiles_wb' 
    base_number = '1120907032.01012'
#    output_iterdb(rhot0,ne,te/e,ni,ti/e,file_out_base+base_number,base_number,'9999',vrot=omega_tor_full,nimp=nz)
#    output_iterdb(rhot0,ne,te/e,ni,tb_full/e,file_out_base+base_number,base_number,'9999',vrot=omega_tor_full,nimp=nz)

    rhop=np.sqrt(psi0)
#    f = open('gene_profiles'+file_out_base+'_gene_e.dat','w')
    f = open('profiles_e.dat','w')
    f.write('# 1.rho_tor 2.rho_pol 3.Te(kev) 4.ne(10^19m^-3)\n#\n')
    np.savetxt(f,np.column_stack((rhot0,rhop,te/e*1E-03,ne*1E-19)))
    f.close()

#    f = open('gene_profiles'+file_out_base+'_gene_i.dat','w')
    f = open('profiles_i.dat','w')
    f.write('# 1.rho_tor 2.rho_pol 3.Ti(kev) 4.ni(10^19m^-3)\n#\n')
#    np.savetxt(f,np.column_stack((rhot0,psi0,ti/e,ni)))
    np.savetxt(f,np.column_stack((rhot0,rhop,tb_full/e*1E-03,ni*1E-19)))
    f.close()
		
    f = open('profiles_z.dat','w')
    f.write('# 1.rho_tor 2.rho_pol 3.Ti(kev) 4.nz(10^19m^-3)\n#\n')
#    np.savetxt(f,np.column_stack((rhot0,psi0,ti/e,ni)))
    np.savetxt(f,np.column_stack((rhot0,rhop,tb_full/e*1E-03,nz*1E-19)))
    f.close()
