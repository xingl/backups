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
#from read_cmod_pfile import *
#from calc_gammaE import *
#from read_cxrs_file import *
#from read_Er import *
from read_iterdb_file import *
#from w_iterdb import *
#from matplotlib.backends.backend_pdf import PdfPages

e = 1.6*10**(-19)
a = 2.6863997038399379E-01

efit_file_name = 'g1120907032.01012'
p_file_name = 'p1120907032.01012'
cxrs_file_name = 'cxrs1120907032.01010_v20140623'
iterdb_filename = 'cmod1120907032.01012.iterdb'

rhot_idb, te_idb, ti_idb, ne_idb, ni_idb, nb_idb, vrot_idb = read_iterdb_file(iterdb_filename)
    
plt.plot(rhot_idb,te_idb,'x',label='te iterdb')
plt.plot(rhot_idb,ti_idb,'x',label='ti iterdb')
plt.xlabel('rhot_n')
plt.legend()
plt.show()

plt.plot(rhot_idb,ne_idb,'x',label='ne iterdb')
plt.plot(rhot_idb,ni_idb,'x',label='ni iterdb')
plt.plot(rhot_idb,nb_idb,'x',label='nb iterdb')
plt.xlabel('rhot_n')
plt.legend()
plt.show()

plt.plot(rhot_idb,vrot_idb,'.',label='iterdb')
plt.axis((0.9,1.0,-150000,20000))
plt.legend()
plt.show()

