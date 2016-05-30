#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import os
import matplotlib.pyplot as plt
from fieldlib import *
from ParIO import * 
import optparse as op
from read_write_geometry import *
from subprocess import call
import sys
import math

parser=op.OptionParser(description='Plots eigenfunction, kperp and omega_d.')
options,args=parser.parse_args()
suffix = args[0]

if suffix == '' :
   suffix = '.dat'
else:
   suffix = '_'+suffix

par = Parameters()
par.Read_Pars('parameters'+suffix)
pars = par.pardict

field = fieldfile('field'+suffix,pars)
field.set_time(field.tfld[-1])

ikx_grid = np.arange(-field.nx/2+1,field.nx/2+1)

theta_grid = np.zeros(field.nx*field.nz,dtype='float128')
phi = np.zeros(field.nx*field.nz,dtype='complex128')

edge_opt = 2.
theta_even = np.linspace(-np.pi,np.pi,field.nz,endpoint=False)
N = np.arcsinh(edge_opt*theta_even[0])/theta_even[0]
theta = 1./edge_opt*np.sinh(N*theta_even)

for i in ikx_grid:
	this_theta_grid = i*2.*np.pi+theta
	theta_grid[(i-ikx_grid[0])*field.nz:(i-ikx_grid[0]+1)*field.nz]=this_theta_grid
delta_z = np.pi*20.
z0 = 18.*np.pi
k0 = np.pi
#phi = np.exp(-(theta_grid-z0)**2/delta_z**2)*np.sin(k0*(theta_grid-z0))+np.exp(-(theta_grid+z0)**2/delta_z**2)*np.sin(k0*(theta_grid+z0))
phi = np.exp(-(theta_grid-z0)**2/delta_z**2)*np.sin(k0*(theta_grid-z0))
plt.plot(theta_grid,np.real(phi),label=r'$Re[\phi]$',color='red')
plt.plot(theta_grid,np.imag(phi),label=r'$Im[\phi]$',color='blue')
plt.plot(theta_grid,np.abs(phi),label=r'$|\phi|$',color='black')
plt.legend()
plt.grid()
plt.show()
zi=complex(0,1)
phi_kz = np.empty(0,dtype='complex128')
kz_grid = np.linspace(-field.nz/4,field.nz/4,1000,endpoint=False)
for k in np.arange(len(kz_grid)):
    print kz_grid[k]
    integral = 0.
    for i in np.arange(len(theta_grid)-1):
        integral = integral + 0.5*(phi[i]*np.exp(zi*kz_grid[k]*theta_grid[i])\
        +phi[i+1]*np.exp(zi*kz_grid[k]*theta_grid[i+1]))*\
        (theta_grid[i+1]-theta_grid[i])
    this_phi_kz = integral/(theta_grid[-1]-theta_grid[0])
    phi_kz = np.append(phi_kz,this_phi_kz)
#kz_grid_prime = kz_grid/pars['major_R']/pars['q0']
#plt.plot(kz_grid_prime,np.abs(phi_kz),label='abs(phi_kz)')
#plt.plot(kz_grid_prime,np.real(phi_kz),label='real(phi_kz)')
#plt.plot(kz_grid_prime,np.imag(phi_kz),label='imag(phi_kz)')
plt.plot(kz_grid,np.abs(phi_kz),label='abs(phi_kz)')
plt.plot(kz_grid,np.real(phi_kz),label='real(phi_kz)')
plt.plot(kz_grid,np.imag(phi_kz),label='imag(phi_kz)')
plt.xlabel('kz')
plt.legend()
plt.show()
