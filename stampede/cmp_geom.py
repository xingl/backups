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
options,args = parser.parse_args()
geom_file_1 = args[0]
geom_file_2 = args[1]
gpars_1,geometry_1 = read_geometry_local(geom_file_1)
gpars_2,geometry_2 = read_geometry_local(geom_file_2)

geom_list = ['gdBdx','gdBdz']
for i in geom_list:
#for i in geometry_1:
    plt.plot(geometry_1[i],label='first geometry')
    plt.plot(geometry_2[i],label='second geometry')
    plt.xlabel('z')
    plt.title(i)
    plt.legend()
    plt.show()
