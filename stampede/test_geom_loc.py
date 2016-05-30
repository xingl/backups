#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from sys import argv,exit,stdout
import optparse as op
from ParIO import *
import os
import sys
import optparse as op
from subprocess import call
from read_write_geometry import *
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

parser = op.OptionParser()
parser.add_option('--write_new','-w',action='store_const',const=1,help='write new tracer_efit file',default=False)
parser.add_option('--cmp_lowres','-l',action='store_const',const=1,help='compare with low res tracer_efit file',default=False)
options,args = parser.parse_args()
#suffix = args[0]
write_new = options.write_new
cmp_lowres = options.cmp_lowres

#geom_file = 'tracer_efit_'+suffix
geom_file = args[0]
lowres_geom_file = args[1]
gene_gpars,gene_geometry = read_geometry_local(geom_file)
if cmp_lowres:
    #lowres_geom_file = 'lowres_tracer_efit_'+suffix
    lowres_gpars,lowres_geometry = read_geometry_local(lowres_geom_file)

geom_list = ['gdBdx','gdBdz']
if (not cmp_lowres) :
    #for i in geometry:
    for i in geom_list:
        info = gene_geometry[i]
        z = np.arange(len(gene_geometry[i]))
        #f=open(i+suffix,'w')
        #np.savetxt(f,np.column_stack((z,info)))
        #f.close()
        plt.plot(gene_geometry[i])
        plt.title(i)
        plt.show()

#new_file_name = 'new_tracer_efit_'+suffix

if write_new:
    gpars,geometry = read_geometry_local(geom_file)

    new_dBdx_file = 'newdBdx'+suffix
    newdBdxfile = np.genfromtxt(new_dBdx_file)
    new_dBdx = newdBdxfile[:,1]
    geometry['gdBdx'] = new_dBdx

    new_dBdz_file = 'newdBdz'+suffix
    newdBdzfile = np.genfromtxt(new_dBdz_file)
    new_dBdz = newdBdzfile[:,1]
    geometry['gdBdz'] = new_dBdz

    write_tracer_efit_file_local(gpars,geometry,new_file_name)

if cmp_lowres:
#    new_gpars,new_geometry = read_geometry_local(new_file_name)
    #for i in geometry:
    for i in geom_list:
        plt.plot(lowres_geometry[i],label='lowres')
        plt.plot(gene_geometry[i],label='highres 901x901')
#        plt.plot(gene_geometry[i],label='highres 257x257')
#        plt.plot(gene_geometry[i],label='highres 513x513')
#        plt.plot(new_geometry[i],label='highres smoothed')
        plt.legend(loc=2)
	#tt='x0 = '+suffix+', '+i
        #plt.title(tt)
        plt.show()

