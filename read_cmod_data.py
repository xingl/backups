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
#psi0, ne, te, ni, ti_pfile, nz_pfile=read_cmod_pfile(p_file_name,shift_Ti=False,shift_Te=False,add_impurity=True)
psipne, ne, psipte, te, psipni, ni, psipti, ti=read_cmod_pfile_raw(p_file_name)



# shift cxrs measurements out by 0.005 psip_n
cxrs_shift = 0.005  #Shift in poloidal flux coord
psip_tb,tb,psip_nb,nb,psip_Er,Er,psip_vpolb,vpolb,psip_vtorb,vtorb = read_cxrs_file_full(cxrs_file_name)
psip_tb = psip_tb+cxrs_shift
psip_nb = psip_nb+cxrs_shift
psip_Er = psip_Er+cxrs_shift
psip_vpolb = psip_vpolb+cxrs_shift
psip_vtorb = psip_vtorb+cxrs_shift

plt.plot(psip_tb,tb,label='shifted tb_cxrs')
plt.plot(psipte,te,label='te_pfile')
plt.plot(psipti,ti,label='ti_pfile')
plt.ylabel('KeV')
plt.legend()
plt.show()

plt.plot(psipne,ne,label='ne_pfile')
plt.plot(psipni,ni,label='ni_pfile')
plt.plot(psip_nb,nb*1E-02,label='shifted nb_cxrs')
plt.ylabel('10^20 m^-3')
plt.legend()
plt.show()

plt.plot(psip_Er,Er,label='Er_cxrs')
plt.ylabel('kV/m')
plt.legend()
plt.show()

plt.plot(psip_vpolb,vpolb,label='vpolb_cxrs')
plt.plot(psip_vtorb,vtorb,label='vtorb_cxrs')
plt.ylabel('km/s')
plt.legend()
plt.show()
