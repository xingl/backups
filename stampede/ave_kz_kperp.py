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

rescale_phi = 2.

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
print 't = ', field.tfld[-1]

if field.nx%2:
   zgrid = np.linspace(-field.nx,field.nx,field.nx*field.nz,endpoint=False)

else:
   zgrid = np.linspace(-(field.nx-1),(field.nx+1),field.nx*field.nz,endpoint=False)

ikx_grid = np.arange(-field.nx/2+1,field.nx/2+1)
#ikx_grid = np.arange(-field.nx/4+1,field.nx/4+1)
kperp = np.zeros(field.nx*field.nz,dtype='float128')
omega_d = np.zeros(field.nx*field.nz,dtype='float128')
theta_grid = np.zeros(field.nx*field.nz,dtype='float128')
phi = np.zeros(field.nx*field.nz,dtype='complex128')
phi2 = np.zeros(field.nx*field.nz,dtype='complex128')

if 'x_local' in pars:
    if pars['x_local']:
        x_local = True
    else:
        x_local = False 
else:
    x_local = True

if 'kx_center' in pars:
    kx_center = pars['kx_center']
else:
    kx_center = 0.

if x_local:
	kymin = pars['kymin']
	lx = pars['lx']
	geom_file = 'tracer_efit'+suffix
	#dkx = 2*np.pi/lx
	dkx = 2*np.pi*pars['shat']*kymin
	gpars,geometry = read_geometry_local(geom_file)

	ggxx = geometry['ggxx']
	ggxy = geometry['ggxy']
	ggyy = geometry['ggyy']
	ggxz = geometry['ggxz']
	ggyz = geometry['ggyz']
	gdBdx = geometry['gdBdx']
	gdBdy = geometry['gdBdy']
	gdBdz = geometry['gdBdz']
	gBfield = geometry['gBfield']

for i in ikx_grid:
	kx = i*dkx+kx_center
	this_kperp = np.sqrt(ggxx*kx**2+2.*ggxy*kx*kymin+ggyy*kymin**2)

	gamma1 = ggxx*ggyy-ggxy**2
	gamma2 = ggxx*ggyz-ggxy*ggxz
	gamma3 = ggxy*ggyz-ggyy*ggxz
	Kx = -(gdBdy+gamma2/gamma1*gdBdz)
	Ky = gdBdx-gamma3/gamma1*gdBdz
	omega_d_y = Ky*kymin
	omega_d_x = Kx*kx
	this_omegad =  (omega_d_y + omega_d_x)/gBfield

	kperp[(i-ikx_grid[0])*field.nz:(i-ikx_grid[0]+1)*field.nz]=this_kperp
	omega_d[(i-ikx_grid[0])*field.nz:(i-ikx_grid[0]+1)*field.nz]=this_omegad

#plt.plot(zgrid,kperp,label='kperp',color='lawngreen')
#plt.plot(zgrid,omega_d,label='omega_d',color='aqua')
#plt.xlabel('z/pi')
#plt.title('kymin = '+str(kymin))
	

if 'n0_global' in pars:
        phase_fac = -np.e**(-2.0*np.pi*(0.0+1.0J)*pars['n0_global']*pars['q0'])
else:
        phase_fac = -1.0
    
if pars['shat'] > 0.:
    for i in ikx_grid:
	this_phi = field.phi()[:,0,i]*phase_fac**i
	phi[(i-ikx_grid[0])*field.nz:(i-ikx_grid[0]+1)*field.nz]=this_phi
else:
    for i in ikx_grid:
	this_phi = field.phi()[:,0,-i]*phase_fac**i
	phi[(i-ikx_grid[0])*field.nz:(i-ikx_grid[0]+1)*field.nz]=this_phi
phi = phi/field.phi()[field.nz/2,0,0]
phi = phi/max(np.abs(phi))*rescale_phi
    

ave_kperp = sum(kperp*abs(phi)**2)/sum(abs(phi)**2)
print 'weighted k_perp', ave_kperp
ave_omegad = sum(omega_d*abs(phi)**2)/sum(abs(phi)**2)
print 'weighted omega_d', ave_omegad

theta_even = np.linspace(-np.pi,np.pi,field.nz,endpoint=False)
#N = np.arcsinh(pars['edge_opt']*theta_even[-1])/theta_even[-1]
N = np.arcsinh(pars['edge_opt']*theta_even[0])/theta_even[0]
theta = 1./pars['edge_opt']*np.sinh(N*theta_even)
#plt.plot(theta_even,'.')
#plt.plot(theta,'.')
#plt.show()

delta_z = 2.6
#phi2[-ikx_grid[0]*field.nz:(-ikx_grid[0]+1)*field.nz] = np.exp(-theta**2/delta_z**2)*rescale_phi

for i in ikx_grid:
#	this_theta_grid = i*2+theta/np.pi
	this_theta_grid = i*2.*np.pi+theta
	theta_grid[(i-ikx_grid[0])*field.nz:(i-ikx_grid[0]+1)*field.nz]=this_theta_grid
#phi2 = np.exp(-theta_grid**2/delta_z**2)*rescale_phi
#plt.plot(theta_grid,'.')
#plt.show()
plt.title(r'$\phi$')
plt.plot(theta_grid,np.real(phi),label=r'$Re[\phi]$',color='red')
plt.plot(theta_grid,np.imag(phi),label=r'$Im[\phi]$',color='blue')
plt.plot(theta_grid,np.abs(phi),label=r'$|\phi|$',color='black')
#plt.plot(theta_grid,np.abs(phi2),label=r'$|\phi|$')
plt.legend()
plt.grid()
##plt.axis([-8.,9.,-2.,5.])
plt.show()

zi=complex(0,1)
phi_kz = np.empty(0,dtype='complex128')
kz_grid = np.linspace(-field.nz/4,field.nz/4,1000,endpoint=False)
#kz_grid = np.linspace(-field.nz/8,field.nz/8,500,endpoint=False)
for k in np.arange(len(kz_grid)):
    #print kz_grid[k]
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
ave_kz_2 = sum(kz_grid**2*abs(phi_kz)**2)/sum(abs(phi_kz)**2)
ave_kz = np.sqrt(ave_kz_2)
print 'weighted kz', ave_kz
