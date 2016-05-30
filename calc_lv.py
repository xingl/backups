#! /usr/bin/python

from pylab import *
from sys import argv,exit,stdout
import optparse as op
import matplotlib.pyplot as plt
import numpy as np
from interp import *

parser = op.OptionParser()
options,args = parser.parse_args()
c_buffer_x = args[0]
buffer_size = args[1]
l_buffer_x = float(c_buffer_x) - 0.4*float(buffer_size)
r_buffer_x = float(c_buffer_x) + 0.4*float(buffer_size)

pdata=np.genfromtxt('p_info.dat')
rhot=pdata[:,0]
te=pdata[:,2]

e = 1.6*10**(-19)

rhot_fine = linspace(rhot[0],rhot[-1],10*len(rhot))
te_fine = interp(rhot,te,rhot_fine)

l_ind = np.argmin(abs(rhot_fine - float(l_buffer_x)))
c_ind = np.argmin(abs(rhot_fine - float(c_buffer_x)))
r_ind = np.argmin(abs(rhot_fine - float(r_buffer_x)))

te_l = te_fine[l_ind]
te_c = te_fine[c_ind]
te_r = te_fine[r_ind]

lv = 3*np.sqrt(te_l/te_c)
lw = lv**2

print 'lv = ', lv
print 'te_l = ', te_l*0.001/e
print 'te_c = ', te_c*0.001/e
print 'te_r = ', te_r*0.001/e
print 'lw = ', lw
print 'buffer', l_buffer_x, r_buffer_x

nv = 48*np.sqrt(te_l/te_r)
nw = 16*te_l/te_r

print 'nv = ', nv
print 'nw = ', nw

plt.plot(rhot,te)
plt.xlabel('rhot_n')
plt.ylabel('te')
plt.axis((0.95,1.0,0,1.3E-16))
plt.axvline(x=l_buffer_x,color='green')
plt.axvline(x=r_buffer_x,color='green')
plt.axvline(x=c_buffer_x,color='red')

plt.show()
