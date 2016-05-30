import numpy as np
import matplotlib.pyplot as plt
from ParIO import *
import os
"""
Looks through scanfiles directories to find gamma_omega_ky.dat and modes_ky.dat then produces x0 vs ky contour plot.
"""

filter_KBMs = False

scan_nums = []

def parse_modes_file(dir):
    f=open(dir,'r')
    temp = f.read()
    lines = temp.split('\n')
    modes = []
    for line in lines:
        if line:
            if line[0] != "#":
                modes.append(line.split()[1])
    print "modes",modes
    return modes

def set_KBMs_to_0(data,modes):
   for x0 in data:
       for i in range(len(data[x0][:,0])): 
           if modes[x0][i] == 'KBM':
               data[x0][i,1]=0.0

   return data



keep_going = True
dirnum = 0

x0_array = np.empty(0)
data = {}
modes = {}
while keep_going:
    if scan_nums:
        if dirnum+1 <= len(scan_nums):
            dnum = int(float(scan_nums[dirnum]))
        else:
            break
    else:
        dnum = dirnum
    str_dirnum = '0000'+str(dnum)
    str_dirnum = str_dirnum[-4:]
    dir = 'scanfiles'+str_dirnum
    if os.path.exists(dir):
        if os.path.exists(dir+'/'+'gamma_omega_ky.dat'):
            f=open(dir+'/'+'gamma_omega_ky.dat','r')
            ftemp = f.read()
            f.close()
            lines = ftemp.split('\n')
            x0 = lines[0].split()[-1]
            x0_array = np.append(x0_array,float(x0))
            data[x0] = np.genfromtxt(dir+'/'+'gamma_omega_ky.dat')
            modes[x0] = parse_modes_file(dir+'/'+'modes_ky.dat')
        dirnum += 1
    else:
        keep_going=False

if filter_KBMs:
    data = set_KBMs_to_0(data,modes)
	              
#Sort x0 array:
x0_array = np.sort(x0_array) 

for i in range(len(x0_array)):
    if i==0:
        data_out = data[str(x0_array[i])][:,1]
    else:
        data_out = np.column_stack((data_out,data[str(x0_array[i])][:,1]))

ky_array = data[str(x0_array[0])][:,0]
#plt.contourf(x0_array,ky_array,data_out,100,vmax=1.5)

plt.contourf(x0_array,ky_array,data_out,100)
plt.colorbar()
for i in modes:
    print "i in modes",i
    for j in range(len(ky_array)):
        if modes[i][j] == 'MTM':
            mark = 'x'
            col = 'red'
        elif modes[i][j] == 'KBM':
            mark = 'o'
            col = 'black'
        elif modes[i][j] == 'ID':
            mark = (5,2,0)
            col = 'magenta'
        elif modes[i][j] == 'NMTM' or modes[i][j] == 'NBP':
            mark = '+'
            col = 'orange'
        else:
            mark = '+'
            col = 'green'
        plt.plot(float(i),ky_array[j],color=col,marker=mark,markeredgewidth=2,markersize=5)
ylim = plt.ylim()
#plt.ylim(ylim[0]-0.025,ylim[1]+0.025)
plt.ylim(ylim[0]-0.025,0.5)
xlim = plt.xlim()
plt.xlim(xlim[0]-0.002,xlim[1]+0.002)
#plt.axes((ax[0]-0.002,ax[1]+0.002,ax[2]-0.025,ax[2]+0.025))
plt.xlabel(r'$\rho_{tor}$',size=18)
plt.ylabel(r'$k_y \rho_s$',size=18)
plt.show()

homedir = os.path.expanduser("~")
thisdir = os.getcwd()
case = thisdir.split('/')[-1]
rbs = np.genfromtxt(homedir+'/pmv_eqs/'+case+'/rbsProfs')

for i in modes:
    print "i in modes",i
    for j in range(len(ky_array)):
        if modes[i][j] == 'MTM':
            mark = 'x'
            col = 'red'
        elif modes[i][j] == 'KBM':
            mark = 'o'
            col = 'black'
        elif modes[i][j] == 'ID':
            mark = (5,2,0)
            col = 'magenta'
        elif modes[i][j] == 'NMTM' or modes[i][j] == 'NBP':
            mark = '+'
            col = 'orange'
        else:
            mark = '+'
            col = 'green'
        if ky_array[j] <= 1.0:
             plt.plot(float(i),data[i][j,1],color=col,marker=mark,markeredgewidth=2,markersize=5)
#plt.ylim(ylim[0]-0.025,ylim[1]+0.025)
plt.plot(rbs[:,0],abs(rbs[:,9])/2.5,label=r'$\gamma_{E \times B}/2.5 (c_s/a)$')
#plt.axes((ax[0]-0.002,ax[1]+0.002,ax[2]-0.025,ax[2]+0.025))
ax=plt.axis()
plt.axis([x0_array[0]-0.02,1.0,0.0,ax[3]])
plt.legend()
plt.xlabel(r'$\rho_{tor}$',size=18)
plt.ylabel(r'$\gamma (c_s/a)$',size=18)
plt.show()





#avg_gr = np.zeros(len(ky_array))
#for i in range(len(ky_array)):
#    num_x0 = 0
#    for j in range(len(x0_array)):
#        if data_out[i,j] == data_out[i,j] and data_out[i,j] != 0.0:
#            avg_gr[i] += data_out[i,j]
#            num_x0 += 1
#    if num_x0 > 0:
#        avg_gr[i] = avg_gr[i] / float(num_x0)    
#    else:
#        avg_gr[i] = 0.0
#
#plt.plot(ky_array,avg_gr,'x-')
#plt.xlabel(r'$k_y \rho_s$',size=18)
#plt.ylabel(r'$\gamma (c_s/a)$',size=18)
#plt.show()
#
#avg_gr_fname = 'gr_x0_avg'
#if filter_KBMs:
#    avg_gr_fname += '_noKBM.dat'
#else:
#    avg_gr_fname += '.dat'
#
#np.savetxt(avg_gr_fname,np.column_stack((ky_array,avg_gr)))




#data_out2 = 1.0*data_out
#print "np.shape(data_out)",np.shape(data_out)
#for i in range(len(ky_array)):
#    for j in range(len(x0_array)):
#        if data_out2[i,j] > 2.0:
#            data_out2[i,j] = 2.0
#        data_out2[i,j] = data_out2[i,j]/ky_array[i]**2
#     
#    
#plt.contourf(x0_array,ky_array,data_out2,100)
#plt.colorbar()
#for i in modes:
#    print "i in modes",i
#    for j in range(len(ky_array)):
#        if modes[i][j] == 'MTM':
#            mark = 'x'
#            col = 'red'
#        elif modes[i][j] == 'KBM':
#            mark = 'o'
#            col = 'black'
#        elif modes[i][j] == 'NMTM' or modes[i][j] == 'NBP':
#            mark = '+'
#            col = 'yellow'
#        else:
#            mark = '+'
#            col = 'green'
#        plt.plot(float(i),ky_array[j],color=col,marker=mark,markeredgewidth=2,markersize=5)
#ylim = plt.ylim()
##plt.ylim(ylim[0]-0.025,ylim[1]+0.025)
#plt.ylim(ylim[0]-0.025,0.5)
#xlim = plt.xlim()
#plt.xlim(xlim[0]-0.002,xlim[1]+0.002)
##plt.axes((ax[0]-0.002,ax[1]+0.002,ax[2]-0.025,ax[2]+0.025))
#plt.xlabel(r'$\rho_{tor}$',size=18)
#plt.ylabel(r'$k_y \rho_s$',size=18)
#plt.show()

