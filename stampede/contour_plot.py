#Takes data from GDT output and plots phi^2 contour plot
import matplotlib.pyplot as plt
import numpy as np
import re
from matplotlib import colors, ticker

parfile='parameters_1'
datfile='contions_1_9.dat'

phisq=np.genfromtxt(datfile)
plt.contourf(np.log10(phisq),200)
plt.show()

f=open(parfile,'r')
par_in=f.read()
#print par_in
par_lines=par_in.split('\n')
for i in range(len(par_lines)):
    test=re.search('^lx',par_lines[i])
    if test:
        line=par_lines[i].split()
        #print line
        lx=float(line[2])
    test=re.search('nx0',par_lines[i])
    if test:
        line=par_lines[i].split()
        nx0=float(line[2])
    test=re.search('kymin',par_lines[i])
    if test:
        line=par_lines[i].split()
        kymin=float(line[2])
    test=re.search('nky0',par_lines[i])
    if test:
        line=par_lines[i].split()
        nky0=float(line[2])

print "lx=",lx
print "nx0=",nx0
print "kymin=",kymin
print "nky0=",nky0

kxmax=(nx0/2-1)*(2*np.pi/lx)
kymax=(nky0-1)*kymin
print "kxmax",kxmax
print "kymax",kymax

dims=np.shape(phisq)
kxgrid=np.arange(dims[0])/(float(dims[0]-1))*kxmax*2-kxmax
print kxgrid
kygrid=np.arange(dims[1])/(float(dims[1]-1))*kymax*2-kymax
print kygrid


plt.contourf(kxgrid,kygrid,(np.transpose(np.log10(phisq))),200)
#plt.contourf(kxgrid,kygrid,np.transpose(phisq),locator=ticker.LogLocator())
plt.xlabel(r'$k_x \rho_s$',size=18)
plt.ylabel(r'$k_y \rho_s$',size=18)
plt.colorbar()
plt.title(r'log($\phi^2$)',size=18)
#plt.title(r'$\phi^2$',size=18)
plt.show()



