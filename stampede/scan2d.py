import re
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from pylab import *
import os
import sys
from get_nrg import get_nrg0
from sys import path
from subprocess import call
GENE_path="/home1/02289/xl3676/gene11/tools/python"
path.append(GENE_path)
script_path="/home1/02289/xl3676/scripts"
path.append(script_path)
from ParIO import Parameters
from calc_gr import calc_gr

plot=True
calc_grs=False
n_cl_args=len(sys.argv)
fix_parfile = True
ql_fluxes = False
if len(sys.argv) > 1:# and str(sys.argv[1])=='noplot':
    for i in range(1,n_cl_args):
        if str(sys.argv[i])=='noplot':
            plot=False
        if str(sys.argv[i])=='calc_grs':
            calc_grs=True

if fix_parfile:
    call(['cp','parameters','par_bak'])
    f=open('parameters','r')
    partemp=f.read()
    f.close()
    partemps=partemp.split('\n')
    for i in range(len(partemps)):
        do_again = True
        while do_again:
            if partemps[i] and partemps[i][-1]==' ':
                partemps[i]=partemps[i][0:-1]
            else:
                do_again = False
        #if not 'scan_dims' in partemps[i]:
        #    partemps[i]=partemps[i][0:-1]
    partemp='\n'.join(partemps)
    f=open('parameters','w')
    f.write(partemp)
    f.close()

par=Parameters()
par.Read_Pars('parameters')
#print "namelists:", par.nmldict
print "pardict:", par.pardict
print par.pardict['scan_dims']
print type(par.pardict['scan_dims'])
par.pardict['scan_dims'] = str(par.pardict['scan_dims'])
scan_dim_split=par.pardict['scan_dims'].split()
print len(scan_dim_split)
nspec=int(float(par.pardict['n_spec']))
print "nspec",nspec
 
spec1=par.pardict['name1'][1:-1]
print "spec1",spec1
if nspec>=2:
    spec2=par.pardict['name2'][1:-1]
    print "spec2",spec2
if nspec==3:
    spec3=par.pardict['name3'][1:-1]
    print "spec3",spec3


if len(scan_dim_split)==1:
    svar1=raw_input("Please enter the scan parameter: ")
    scan_dims=np.empty(1)
    scan_dims[0]=int(float(par.pardict['scan_dims'].split()[0]))
    print "Dimensions of scan:",'[',svar1,']=',scan_dims
    
    #svar1=input("Please enter first scan parameter: ")
    #svar2=input("Please enter second scan parameter: ")

    file_num=1
    grid1=np.empty(0)
    gamma=np.empty(0)
    omega=np.empty(0)
    Q1ES=np.empty(0)
    Q2ES=np.empty(0)
    Q3ES=np.empty(0)
    Q1EM=np.empty(0)
    Q2EM=np.empty(0)
    Q3EM=np.empty(0)
    G1ES=np.empty(0)
    G2ES=np.empty(0)
    G3ES=np.empty(0)
    G1EM=np.empty(0)
    G2EM=np.empty(0)
    G3EM=np.empty(0)
    #par=read_parameters('')

    #while still_files:
    ntot=int(scan_dims[0])
    for i in range(ntot):
        file_num=i+1
        suffix=str(file_num)
        suffix='000'+suffix
        suffix=suffix[-4:]
        suffix='_'+suffix
        if not os.path.isfile('parameters'+suffix):
            #still_files=False
            print "file does not exist:",'parameters'+suffix
            print "Inserting nan's"
            grid1=np.append(grid1,float('nan'))
            print "grid1",i
            gamma=np.append(gamma,float('nan'))
            omega=np.append(omega,float('nan'))
            if nspec>=1 and ql_fluxes:
                Q1ES=np.append(Q1ES,float('nan'))
                Q1EM=np.append(Q1EM,float('nan'))
                G1ES=np.append(G1ES,float('nan'))
                G1EM=np.append(G1EM,float('nan'))
            elif nspec>=2 and ql_fluxes:
                Q2ES=np.append(Q2ES,float('nan'))
                Q2EM=np.append(Q2EM,float('nan'))
                G2ES=np.append(G2ES,float('nan'))
                G2EM=np.append(G2EM,float('nan'))
            elif nspec==3 and ql_fluxes:
                Q3ES=np.append(Q3ES,float('nan'))
                Q3EM=np.append(Q3EM,float('nan'))
                G3ES=np.append(G3ES,float('nan'))
                G3EM=np.append(G3EM,float('nan'))

            elif ql_fluxes:
                print "not ready for nspec>2!!!"
                stop
        else:
            par=Parameters()
            par.Read_Pars('parameters'+suffix)
            #scan_dims=par.pardict['scan_dims']
            grid1=np.append(grid1,par.pardict[svar1])
            if os.path.isfile('omega'+suffix):
                print suffix
                frequency=np.genfromtxt('omega'+suffix)    
                if(frequency.any() and frequency[1] != 0.0):
                    gamma=np.append(gamma,frequency[1])
                    omega=np.append(omega,frequency[2])
                    if nspec==1 and ql_fluxes:
                        tn,nrg1=get_nrg0(suffix,nspec=nspec)
                        Q1ES=np.append(Q1ES,nrg1[-1,6]/nrg1[-1,0])
                        Q1EM=np.append(Q1EM,nrg1[-1,7]/nrg1[-1,0])
                        G1ES=np.append(G1ES,nrg1[-1,4]/nrg1[-1,0])
                        G1EM=np.append(G1EM,nrg1[-1,5]/nrg1[-1,0])
                    elif nspec==2 and ql_fluxes:
                        tn,nrg1,nrg2=get_nrg0(suffix,nspec=nspec)
                        Q1ES=np.append(Q1ES,nrg1[-1,6]/nrg1[-1,0])
                        Q1EM=np.append(Q1EM,nrg1[-1,7]/nrg1[-1,0])
                        G1ES=np.append(G1ES,nrg1[-1,4]/nrg1[-1,0])
                        G1EM=np.append(G1EM,nrg1[-1,5]/nrg1[-1,0])
                        Q2ES=np.append(Q2ES,nrg2[-1,6]/nrg2[-1,0])
                        Q2EM=np.append(Q2EM,nrg2[-1,7]/nrg2[-1,0])
                        G2ES=np.append(G2ES,nrg2[-1,4]/nrg2[-1,0])
                        G2EM=np.append(G2EM,nrg2[-1,5]/nrg2[-1,0])
                    elif nspec==3 and ql_fluxes:
                        tn,nrg1,nrg2,nrg3=get_nrg0(suffix,nspec=nspec)
                        Q1ES=np.append(Q1ES,nrg1[-1,6]/nrg1[-1,0])
                        Q1EM=np.append(Q1EM,nrg1[-1,7]/nrg1[-1,0])
                        G1ES=np.append(G1ES,nrg1[-1,4]/nrg1[-1,0])
                        G1EM=np.append(G1EM,nrg1[-1,5]/nrg1[-1,0])
                        Q2ES=np.append(Q2ES,nrg2[-1,6]/nrg2[-1,0])
                        Q2EM=np.append(Q2EM,nrg2[-1,7]/nrg2[-1,0])
                        G2ES=np.append(G2ES,nrg2[-1,4]/nrg2[-1,0])
                        G2EM=np.append(G2EM,nrg2[-1,5]/nrg2[-1,0])
                        Q3ES=np.append(Q3ES,nrg3[-1,6]/nrg3[-1,0])
                        Q3EM=np.append(Q3EM,nrg3[-1,7]/nrg3[-1,0])
                        G3ES=np.append(G3ES,nrg3[-1,4]/nrg3[-1,0])
                        G3EM=np.append(G3EM,nrg3[-1,5]/nrg3[-1,0])
                    elif ql_fluxes:
                        print "not ready for nspec>3!!!"
                        stop

                else:
 		    if calc_grs:
                        gr_temp=calc_gr(suffix,nspec=nspec)
                    else:
                        gr_temp=-1
                    if gr_temp==-1:
                        gamma=np.append(gamma,float('nan'))
                    else:
                        gamma=np.append(gamma,gr_temp)
                    omega=np.append(omega,float('nan'))
                    if nspec==1 and ql_fluxes:
                        Q1ES=np.append(Q1ES,float('nan'))
                        Q1EM=np.append(Q1EM,float('nan'))
                        G1ES=np.append(G1ES,float('nan'))
                        G1EM=np.append(G1EM,float('nan'))
                    elif nspec==2 and ql_fluxes:
                        Q1ES=np.append(Q1ES,float('nan'))
                        Q1EM=np.append(Q1EM,float('nan'))
                        G1ES=np.append(G1ES,float('nan'))
                        G1EM=np.append(G1EM,float('nan'))
                        Q2ES=np.append(Q2ES,float('nan'))
                        Q2EM=np.append(Q2EM,float('nan'))
                        G2ES=np.append(G2ES,float('nan'))
                        G2EM=np.append(G2EM,float('nan'))
                    elif ql_fluxes:
                        print "not ready for nspec>2!!!"
                        stop


            else:
                print suffix, " nan"
                gamma=np.append(gamma,float('nan'))
                omega=np.append(omega,float('nan'))
                print suffix, " nan"
                gamma=np.append(gamma,float('nan'))
                omega=np.append(omega,float('nan'))
                if nspec>=1 and ql_fluxes:
                    Q1ES=np.append(Q1ES,float('nan'))
                    Q1EM=np.append(Q1EM,float('nan'))
                    G1ES=np.append(G1ES,float('nan'))
                    G1EM=np.append(G1EM,float('nan'))
                elif nspec>=2 and ql_fluxes:
                    Q2ES=np.append(Q2ES,float('nan'))
                    Q2EM=np.append(Q2EM,float('nan'))
                    G2ES=np.append(G2ES,float('nan'))
                    G2EM=np.append(G2EM,float('nan'))
                elif nspec==3 and ql_fluxes:
                    Q3ES=np.append(Q3ES,float('nan'))
                    Q3EM=np.append(Q3EM,float('nan'))
                    G3ES=np.append(G3ES,float('nan'))
                    G3EM=np.append(G3EM,float('nan'))
                elif ql_fluxes:
                    print "not ready for nspec>2!!!"
                    stop



    f=open('omega_'+svar1,'w')
    #f.write('#'+svar1+', '+'gamma, omega')
    #f.write('\n')
    print grid1
    print gamma
    if plot:
        plt.plot(grid1,gamma,'x',markeredgewidth=2,markersize=6)
        plt.xlabel(svar1,size=18)
        plt.ylabel(r'$\gamma (c_s/a)$',size=18)
        plt.show()
    temp=np.empty([3,scan_dims[0]])
    temp[0,:]=grid1
    temp[1,:]=gamma
    temp[2,:]=omega
    #print "temp",temp
    np.savetxt('omega_'+svar1,np.transpose(temp))
    f.close()

    if ql_fluxes:
        f=open('qlflux_'+spec1+'_'+svar1,'w')
        f.write('#'+svar1+", GES, GEM, QES, QEM \n")
        temp=np.empty([5,scan_dims[0]])
        temp[0,:]=grid1
        temp[1,:]=G1ES
        temp[2,:]=G1EM
        temp[3,:]=Q1ES
        temp[4,:]=Q1EM
        np.savetxt(f,np.transpose(temp))
        f.close()
        if nspec >= 2:
            f=open('qlflux_'+spec2+'_'+svar1,'w')
            f.write('#'+svar1+", GES, GEM, QES, QEM \n")
            temp=np.empty([5,scan_dims[0]])
            temp[0,:]=grid1
            temp[1,:]=G2ES
            temp[2,:]=G2EM
            temp[3,:]=Q2ES
            temp[4,:]=Q2EM
            np.savetxt(f,np.transpose(temp))
            f.close()
        if nspec >= 3:
            f=open('qlflux_'+spec3+'_'+svar1,'w')
            f.write('#'+svar1+", GES, GEM, QES, QEM \n")
            temp=np.empty([5,scan_dims[0]])
            temp[0,:]=grid1
            temp[1,:]=G3ES
            temp[2,:]=G3EM
            temp[3,:]=Q3ES
            temp[4,:]=Q3EM
            np.savetxt(f,np.transpose(temp))
            f.close()


    if plot:
        plt.plot(grid1,omega,'x',markeredgewidth=2,markersize=6)
        plt.xlabel(svar1,size=18)
        plt.ylabel(r'$\omega (c_s/a)$',size=18)
        plt.show()

    if plot and ql_fluxes:
        plt.plot(grid1,Q1ES,'x',markeredgewidth=2,markersize=6,label='Q'+'('+spec1+')'+'ES')
        plt.plot(grid1,Q1EM,'x',markeredgewidth=2,markersize=6,label='Q'+'('+spec1+')'+'EM')
        plt.plot(grid1,G1ES,'x',markeredgewidth=2,markersize=6,label='G'+'('+spec1+')'+'ES')
        plt.plot(grid1,G1EM,'x',markeredgewidth=2,markersize=6,label='G'+'('+spec1+')'+'EM')
        plt.xlabel(svar1,size=18)
        plt.legend()
        plt.show()
        if nspec>=2:
            plt.plot(grid1,Q2ES,'x',markeredgewidth=2,markersize=6,label='Q'+'('+spec2+')'+'ES')
            plt.plot(grid1,Q2EM,'x',markeredgewidth=2,markersize=6,label='Q'+'('+spec2+')'+'EM')
            plt.plot(grid1,G2ES,'x',markeredgewidth=2,markersize=6,label='G'+'('+spec2+')'+'ES')
            plt.plot(grid1,G2EM,'x',markeredgewidth=2,markersize=6,label='G'+'('+spec2+')'+'EM')
            plt.xlabel(svar1,size=18)
            plt.legend()
            plt.show()
        if nspec>=3:
            plt.plot(grid1,Q3ES,'x',markeredgewidth=2,markersize=6,label='Q'+'('+spec3+')'+'ES')
            plt.plot(grid1,Q3EM,'x',markeredgewidth=2,markersize=6,label='Q'+'('+spec3+')'+'EM')
            plt.plot(grid1,G3ES,'x',markeredgewidth=2,markersize=6,label='G'+'('+spec3+')'+'ES')
            plt.plot(grid1,G3EM,'x',markeredgewidth=2,markersize=6,label='G'+'('+spec3+')'+'EM')
            plt.xlabel(svar1,size=18)
            plt.legend()
            plt.show()

    #if plot:
    #  for i in range(scan_dims[1]):
    #    #print grid1
    #    #print gamma[:,i]
    #    plt.plot(grid1,gamma[:,i],'x',markeredgewidth=2,markersize=6,label=svar2+'='+str(grid2[i]))

    #  plt.xlabel(svar1,size=18)
    #  plt.ylabel(r'$\gamma (c_s/a)$',size=18)
    #  plt.legend(loc='upper left')
    #  plt.show()

    #if plot:
    #  for i in range(scan_dims[1]):
    #    #print grid1
    #    #print omega[:,i]
    #    plt.plot(grid1,omega[:,i],'x',markeredgewidth=2,markersize=6,label=svar2+'='+str(grid2[i]))

    #  plt.xlabel(svar1,size=18)
    #  plt.ylabel(r'$\omega (c_s/a)$',size=18)
    #  plt.legend(loc='lower left')
    #  plt.show()


elif len(scan_dim_split)==2:

    #svar1=input("Please enter first scan parameter: ")
    #svar2=input("Please enter second scan parameter: ")
    print par.pardict
    svar2=raw_input("Please enter first scan parameter: ")
    svar1=raw_input("Please enter second scan parameter: ")
    scan_dims=np.empty(2)
    scan_dims[0]=int(float(par.pardict['scan_dims'].split()[0]))
    scan_dims[1]=int(float(par.pardict['scan_dims'].split()[1]))
    print "Dimensions of scan:",'[',svar1,',',svar2,']=',scan_dims
    

    still_files=True
    file_num=1
    grid1=np.empty(0)
    grid2=np.empty(0)
    gamma=np.empty(0,dtype='float64')
    omega=np.empty(0,dtype='float64')
    Q1ES=np.empty(0)
    Q2ES=np.empty(0)
    Q3ES=np.empty(0)
    Q1EM=np.empty(0)
    Q2EM=np.empty(0)
    Q3EM=np.empty(0)
    G1ES=np.empty(0)
    G2ES=np.empty(0)
    G3ES=np.empty(0)
    G1EM=np.empty(0)
    G2EM=np.empty(0)
    G3EM=np.empty(0)

    #par=read_parameters('')

    #while still_files:
    ntot=int(scan_dims[0]*scan_dims[1])
    for i in range(ntot):
        file_num=i+1
        suffix=str(file_num)
        suffix='000'+suffix
        suffix=suffix[-4:]
        suffix='_'+suffix
        if not os.path.isfile('parameters'+suffix):
            #still_files=False
            print "file does not exist:",'parameters'+suffix
            print "Inserting nan's"
            grid1=np.append(grid1,float('nan'))
            grid2=np.append(grid2,float('nan'))
            gamma=np.append(gamma,float('nan'))
            omega=np.append(omega,float('nan'))
            if nspec>=1 and ql_fluxes:
                Q1ES=np.append(Q1ES,float('nan'))
                Q1EM=np.append(Q1EM,float('nan'))
                G1ES=np.append(G1ES,float('nan'))
                G1EM=np.append(G1EM,float('nan'))
            elif nspec>=2 and ql_fluxes:
                Q2ES=np.append(Q2ES,float('nan'))
                Q2EM=np.append(Q2EM,float('nan'))
                G2ES=np.append(G2ES,float('nan'))
                G2EM=np.append(G2EM,float('nan'))
            elif nspec==3 and ql_fluxes:
                Q3ES=np.append(Q3ES,float('nan'))
                Q3EM=np.append(Q3EM,float('nan'))
                G3ES=np.append(G3ES,float('nan'))
                G3EM=np.append(G3EM,float('nan'))
            elif ql_fluxes:
                print "not ready for nspec>2!!!"
                stop

        else:
            par=Parameters()
            par.Read_Pars('parameters'+suffix)
            #print par.pardict
            #scan_dims=par.pardict['scan_dims']
            grid1=np.append(grid1,par.pardict[svar1])
            grid2=np.append(grid2,par.pardict[svar2])
            if os.path.isfile('omega'+suffix):
                print suffix
                frequency=np.genfromtxt('omega'+suffix)    
                if(frequency.any() and frequency[1] != 0.0):
                    print suffix, frequency[1]
                    gamma=np.append(gamma,float(frequency[1]))
                    omega=np.append(omega,float(frequency[2]))
                    print suffix, gamma[-1]
                    if nspec==1 and ql_fluxes:
                        tn,nrg1=get_nrg0(suffix,nspec=nspec)
                        Q1ES=np.append(Q1ES,nrg1[-1,6]/nrg1[-1,0])
                        Q1EM=np.append(Q1EM,nrg1[-1,7]/nrg1[-1,0])
                        G1ES=np.append(G1ES,nrg1[-1,4]/nrg1[-1,0])
                        G1EM=np.append(G1EM,nrg1[-1,5]/nrg1[-1,0])
                    elif nspec==2 and ql_fluxes:
                        tn,nrg1,nrg2=get_nrg0(suffix,nspec=nspec)
                        Q1ES=np.append(Q1ES,nrg1[-1,6]/nrg1[-1,0])
                        Q1EM=np.append(Q1EM,nrg1[-1,7]/nrg1[-1,0])
                        G1ES=np.append(G1ES,nrg1[-1,4]/nrg1[-1,0])
                        G1EM=np.append(G1EM,nrg1[-1,5]/nrg1[-1,0])
                        Q2ES=np.append(Q2ES,nrg2[-1,6]/nrg2[-1,0])
                        Q2EM=np.append(Q2EM,nrg2[-1,7]/nrg2[-1,0])
                        G2ES=np.append(G2ES,nrg2[-1,4]/nrg2[-1,0])
                        G2EM=np.append(G2EM,nrg2[-1,5]/nrg2[-1,0])
                    elif nspec==3 and ql_fluxes:
                        tn,nrg1,nrg2,nrg3=get_nrg0(suffix,nspec=nspec)
                        Q1ES=np.append(Q1ES,nrg1[-1,6]/nrg1[-1,0])
                        Q1EM=np.append(Q1EM,nrg1[-1,7]/nrg1[-1,0])
                        G1ES=np.append(G1ES,nrg1[-1,4]/nrg1[-1,0])
                        G1EM=np.append(G1EM,nrg1[-1,5]/nrg1[-1,0])
                        Q2ES=np.append(Q2ES,nrg2[-1,6]/nrg2[-1,0])
                        Q2EM=np.append(Q2EM,nrg2[-1,7]/nrg2[-1,0])
                        G2ES=np.append(G2ES,nrg2[-1,4]/nrg2[-1,0])
                        G2EM=np.append(G2EM,nrg2[-1,5]/nrg2[-1,0])
                        Q3ES=np.append(Q3ES,nrg3[-1,6]/nrg3[-1,0])
                        Q3EM=np.append(Q3EM,nrg3[-1,7]/nrg3[-1,0])
                        G3ES=np.append(G3ES,nrg3[-1,4]/nrg3[-1,0])
                        G3EM=np.append(G3EM,nrg3[-1,5]/nrg3[-1,0])
                    elif ql_fluxes:
                        print "not ready for nspec>3!!!"
                        stop

                else:
                    if calc_grs:
                        print "omega file is empty.  Attempting to calculate growth rate from nrg file."
                        gr_temp=calc_gr(suffix,nspec=nspec)
                    else:
                        gr_temp=-1
                    if gr_temp != -1:
                        gamma=np.append(gamma,gr_temp)
                    else:
                        gamma=np.append(gamma,float('nan'))
                    omega=np.append(omega,float('nan'))
                    if nspec==1 and ql_fluxes:
                        Q1ES=np.append(Q1ES,float('nan'))
                        Q1EM=np.append(Q1EM,float('nan'))
                        G1ES=np.append(G1ES,float('nan'))
                        G1EM=np.append(G1EM,float('nan'))
                    elif nspec==2 and ql_fluxes:
                        Q1ES=np.append(Q1ES,float('nan'))
                        Q1EM=np.append(Q1EM,float('nan'))
                        G1ES=np.append(G1ES,float('nan'))
                        G1EM=np.append(G1EM,float('nan'))
                        Q2ES=np.append(Q2ES,float('nan'))
                        Q2EM=np.append(Q2EM,float('nan'))
                        G2ES=np.append(G2ES,float('nan'))
                        G2EM=np.append(G2EM,float('nan'))
                    elif nspec==3 and ql_fluxes:
                        Q1ES=np.append(Q1ES,float('nan'))
                        Q1EM=np.append(Q1EM,float('nan'))
                        G1ES=np.append(G1ES,float('nan'))
                        G1EM=np.append(G1EM,float('nan'))
                        Q2ES=np.append(Q2ES,float('nan'))
                        Q2EM=np.append(Q2EM,float('nan'))
                        G2ES=np.append(G2ES,float('nan'))
                        G2EM=np.append(G2EM,float('nan'))
                        Q3ES=np.append(Q3ES,float('nan'))
                        Q3EM=np.append(Q3EM,float('nan'))
                        G3ES=np.append(G3ES,float('nan'))
                        G3EM=np.append(G3EM,float('nan'))
                    elif ql_fluxes:
                        print "not ready for nspec>2!!!"
                        stop

            else:
                print suffix, " nan"
                gamma=np.append(gamma,float('nan'))
                omega=np.append(omega,float('nan'))
                if nspec==1 and ql_fluxes:
                    Q1ES=np.append(Q1ES,float('nan'))
                    Q1EM=np.append(Q1EM,float('nan'))
                    G1ES=np.append(G1ES,float('nan'))
                    G1EM=np.append(G1EM,float('nan'))
                elif nspec==2 and ql_fluxes:
                    Q1ES=np.append(Q1ES,float('nan'))
                    Q1EM=np.append(Q1EM,float('nan'))
                    G1ES=np.append(G1ES,float('nan'))
                    G1EM=np.append(G1EM,float('nan'))
                    Q2ES=np.append(Q2ES,float('nan'))
                    Q2EM=np.append(Q2EM,float('nan'))
                    G2ES=np.append(G2ES,float('nan'))
                    G2EM=np.append(G2EM,float('nan'))
                elif nspec==3 and ql_fluxes:
                    Q1ES=np.append(Q1ES,float('nan'))
                    Q1EM=np.append(Q1EM,float('nan'))
                    G1ES=np.append(G1ES,float('nan'))
                    G1EM=np.append(G1EM,float('nan'))
                    Q2ES=np.append(Q2ES,float('nan'))
                    Q2EM=np.append(Q2EM,float('nan'))
                    G2ES=np.append(G2ES,float('nan'))
                    G2EM=np.append(G2EM,float('nan'))
                    Q3ES=np.append(Q3ES,float('nan'))
                    Q3EM=np.append(Q3EM,float('nan'))
                    G3ES=np.append(G3ES,float('nan'))
                    G3EM=np.append(G3EM,float('nan'))
                elif ql_fluxes:
                    print "not ready for nspec>2!!!"
                    stop


    #if len(grid1) != ntot:
    #    print "Scan not finished:",len(grid1),"of",scan_dims[0]*scan_dims[1]
    #    while file_num-1 <= ntot:
    #print svar1,grid1
    #print svar2,grid2
    #print "gamma",gamma
    #print "omega",omega
    print "grid1",grid1
    print "len(grid1),len(grid2),len(gamma),len(omega)"
    print len(grid1),len(grid2),len(gamma),len(omega)
    print 'scan_dims', scan_dims
    grid1=np.reshape(grid1,scan_dims)   
    grid1=grid1[:,0]
    #print grid1
    grid2=np.reshape(grid2,scan_dims)   
    grid2=grid2[0,:]
    #print grid2
    gamma=np.reshape(gamma,scan_dims)   
    omega=np.reshape(omega,scan_dims)   
    if ql_fluxes:
        Q1ES=np.reshape(Q1ES,scan_dims)
        Q1EM=np.reshape(Q1EM,scan_dims)
        G1ES=np.reshape(G1ES,scan_dims)
        G1EM=np.reshape(G1EM,scan_dims)
        if nspec>=2:
            Q2ES=np.reshape(Q2ES,scan_dims)
            Q2EM=np.reshape(Q2EM,scan_dims)
            G2ES=np.reshape(G2ES,scan_dims)
            G2EM=np.reshape(G2EM,scan_dims)
        if nspec==3:
            Q3ES=np.reshape(Q3ES,scan_dims)
            Q3EM=np.reshape(Q3EM,scan_dims)
            G3ES=np.reshape(G3ES,scan_dims)
            G3EM=np.reshape(G3EM,scan_dims)

    print "grid1",grid1
    for i in range(int(scan_dims[1])):
        f=open('omega_'+svar2+'_'+str(i),'w')
        #f.write('#'+svar1+', '+'gamma, omega')
        #f.write('\n')
        #print grid1
        #print gamma[:,i]
        if plot:
            plt.title(svar2+'='+str(grid2[i]))
            plt.plot(grid1,gamma[:,i],'x',markeredgewidth=2,markersize=6)
            plt.xlabel(svar1,size=18)
            plt.ylabel(r'$\gamma (c_s/a)$',size=18)
            plt.show()
        temp=np.empty([3,scan_dims[0]])
        temp[0,:]=grid1
        temp[1,:]=gamma[:,i]
        temp[2,:]=omega[:,i]
        print "temp",temp
        np.savetxt('omega_'+svar2+'_'+str(i),np.transpose(temp))
        f.close()

    for i in range(int(scan_dims[1])):
        if ql_fluxes:
            f=open('qlflux_'+spec1+'_'+svar2+'_'+str(i),'w')
            f.write('#'+svar1+", GES, GEM, QES, QEM \n")
            temp=np.empty([5,scan_dims[0]])
            temp[0,:]=grid1
            temp[1,:]=G1ES[:,i]
            temp[2,:]=G1EM[:,i]
            temp[3,:]=Q1ES[:,i]
            temp[4,:]=Q1EM[:,i]
            np.savetxt(f,np.transpose(temp))
            f.close()
            if nspec >=2:
                f=open('qlflux_'+spec2+'_'+svar2+'_'+str(i),'w')
                f.write('#'+svar1+", GES, GEM, QES, QEM \n")
                temp=np.empty([5,scan_dims[0]])
                temp[0,:]=grid1
                temp[1,:]=G2ES[:,i]
                temp[2,:]=G2EM[:,i]
                temp[3,:]=Q2ES[:,i]
                temp[4,:]=Q2EM[:,i]
                np.savetxt(f,np.transpose(temp))
                f.close()
            if nspec ==3:
                f=open('qlflux_'+spec3+'_'+svar2+'_'+str(i),'w')
                f.write('#'+svar1+", GES, GEM, QES, QEM \n")
                temp=np.empty([5,scan_dims[0]])
                temp[0,:]=grid1
                temp[1,:]=G3ES[:,i]
                temp[2,:]=G3EM[:,i]
                temp[3,:]=Q3ES[:,i]
                temp[4,:]=Q3EM[:,i]
                np.savetxt(f,np.transpose(temp))
                f.close()


            f=open('omega_'+svar2+'_'+str(i),'w')
            #f.write('#'+svar1+', '+'gamma, omega')
            #f.write('\n')
            #print grid1
            #print gamma[:,i]
            if plot:
                plt.title(svar2+'='+str(grid2[i]))
                plt.plot(grid1,gamma[:,i],'x',markeredgewidth=2,markersize=6)
                plt.xlabel(svar1,size=18)
                plt.ylabel(r'$\gamma (c_s/a)$',size=18)
                plt.show()
            temp=np.empty([3,scan_dims[0]])
            temp[0,:]=grid1
            temp[1,:]=gamma[:,i]
            temp[2,:]=omega[:,i]
            print "temp",temp
            np.savetxt('omega_'+svar2+'_'+str(i),np.transpose(temp))
            f.close()


    if plot:
      for i in range(int(scan_dims[1])):
        #print grid1
        #print omega[:,i]
        plt.title(svar2+'='+str(grid2[i]))
        plt.plot(grid1,omega[:,i],'x',markeredgewidth=2,markersize=6)
        plt.xlabel(svar1,size=18)
        plt.ylabel(r'$\omega (c_s/a)$',size=18)
        plt.show()

    if plot:
      for i in range(int(scan_dims[1])):
        #print grid1
        #print gamma[:,i]
        plt.plot(grid1,gamma[:,i],'x',markeredgewidth=2,markersize=6,label=svar2+'='+str(grid2[i]))

      plt.xlabel(svar1,size=18)
      plt.ylabel(r'$\gamma (c_s/a)$',size=18)
      plt.legend(loc='upper right')
      plt.show()

    if plot:
      for i in range(int(scan_dims[1])):
        #print grid1
        #print omega[:,i]
        plt.plot(grid1,omega[:,i],'x',markeredgewidth=2,markersize=6,label=svar2+'='+str(grid2[i]))

      plt.xlabel(svar1,size=18)
      plt.ylabel(r'$\omega (c_s/a)$',size=18)
      plt.legend(loc='lower right')
      plt.show()

    if plot and ql_fluxes:
      for i in range(int(scan_dims[1])):
        plt.title(svar2+'='+str(grid2[i]))
        plt.plot(grid1,Q1ES[:,i],'x',markeredgewidth=2,markersize=6,label='Q'+'('+spec1+')'+'ES')
        plt.plot(grid1,Q1EM[:,i],'x',markeredgewidth=2,markersize=6,label='Q'+'('+spec1+')'+'EM')
        plt.plot(grid1,G1ES[:,i],'x',markeredgewidth=2,markersize=6,label='G'+'('+spec1+')'+'ES')
        plt.plot(grid1,G1EM[:,i],'x',markeredgewidth=2,markersize=6,label='G'+'('+spec1+')'+'EM')
        plt.xlabel(svar1,size=18)
        plt.legend()
        plt.show()
        if nspec>=2:
            plt.title(svar2+'='+str(grid2[i]))
            plt.plot(grid1,Q2ES[:,i],'x',markeredgewidth=2,markersize=6,label='Q'+'('+spec2+')'+'ES')
            plt.plot(grid1,Q2EM[:,i],'x',markeredgewidth=2,markersize=6,label='Q'+'('+spec2+')'+'EM')
            plt.plot(grid1,G2ES[:,i],'x',markeredgewidth=2,markersize=6,label='G'+'('+spec2+')'+'ES')
            plt.plot(grid1,G2EM[:,i],'x',markeredgewidth=2,markersize=6,label='G'+'('+spec2+')'+'EM')
            plt.xlabel(svar1,size=18)
            plt.legend()
            plt.show()
        if nspec==3:
            plt.title(svar2+'='+str(grid2[i]))
            plt.plot(grid1,Q3ES[:,i],'x',markeredgewidth=2,markersize=6,label='Q'+'('+spec3+')'+'ES')
            plt.plot(grid1,Q3EM[:,i],'x',markeredgewidth=2,markersize=6,label='Q'+'('+spec3+')'+'EM')
            plt.plot(grid1,G3ES[:,i],'x',markeredgewidth=2,markersize=6,label='G'+'('+spec3+')'+'ES')
            plt.plot(grid1,G3EM[:,i],'x',markeredgewidth=2,markersize=6,label='G'+'('+spec3+')'+'EM')
            plt.xlabel(svar1,size=18)
            plt.legend()
            plt.show()





