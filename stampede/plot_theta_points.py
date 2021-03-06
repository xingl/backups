#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from read_write_geometry import *

geomfile = 'tracer_efit_nITER_b80_3a_98_eo2'


parameters, geom = read_geometry_local(geomfile)
R = geom['gl_R']
Z = geom['gl_z']
B = geom['gBfield']*parameters['Bref']

for i in range(len(R)):
    plt.scatter(R[i],Z[i])
plt.xlabel('R(m)')
plt.ylabel('Z(m)')
plt.show()


