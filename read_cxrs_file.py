import matplotlib.pyplot as plt
import numpy as np
import re

def read_cxrs_file(cxrs_filename):

    cxrs_file=open(cxrs_filename,'r')

    data_in=cxrs_file.read()
    data_linesplit=data_in.split('\n')

    keep_going=1
    i=0
    while keep_going:
        test=re.search('psinorm',data_linesplit[i])
        if test:
            quantity=data_linesplit[i].split()[1]
            print "Reading :",quantity
            if quantity=='T_{Boron5+}(keV)':
                psip_tb=np.empty(0)
                tb=np.empty(0)
                tb_error=np.empty(0)
	        r = 1
		l = 2
                while r:
                    str_temp=data_linesplit[i+l].split()
                    if not len(data_linesplit[i+l].split())==0:
                        psip_tb=np.append(psip_tb,float(str_temp[0]))
                        tb=np.append(tb,float(str_temp[1]))
                        tb_error=np.append(tb_error,float(str_temp[2]))
                        l=l+1
                        i=i+1
                    else:
                        r= (1==2)
            if quantity=='n_{Boron5+}(10^18':
                psip_nb=np.empty(0)
                nb=np.empty(0)
                nb_error=np.empty(0)
                r = 1
                l = 2
                while r:
                    str_temp=data_linesplit[i+l].split()
                    if not len(data_linesplit[i+l].split())==0:
                        psip_nb=np.append(psip_nb,float(str_temp[0]))
                        nb=np.append(nb,float(str_temp[1]))
                        nb_error=np.append(nb_error,float(str_temp[2]))
                        l=l+1
                        i=i+1
                    else:
                        r= (1==2)
            if quantity=='vpol_{Boron5+}(km/s)':
                psip_vpolb=np.empty(0)
                vpolb=np.empty(0)
                vpolb_error=np.empty(0)
                r = 1
                l = 2
                while r:
                    str_temp=data_linesplit[i+l].split()
                    if not len(data_linesplit[i+l].split())==0:
                        psip_vpolb=np.append(psip_vpolb,float(str_temp[0]))
                        vpolb=np.append(vpolb,float(str_temp[1]))
                        vpolb_error=np.append(vpolb_error,float(str_temp[2]))
                        l=l+1
                        i=i+1
                    else:
                        r= (1==2)
            if quantity=='vtor_{Boron5+}(km/s)':
                psip_vtorb=np.empty(0)
                vtorb=np.empty(0)
                vtorb_error=np.empty(0)
                r = 1
                l = 2
                while r:
                    str_temp=data_linesplit[i+l].split()
                    if not len(data_linesplit[i+l].split())==0:
                        psip_vtorb=np.append(psip_vtorb,float(str_temp[0]))
                        vtorb=np.append(vtorb,float(str_temp[1]))
                        vtorb_error=np.append(vtorb_error,float(str_temp[2]))
                        l=l+1
                        i=i+1
                    else:
                        r= (1==2)
            if quantity=='Er(kV/m)':
                psip_Er=np.empty(0)
                Er=np.empty(0)
                Er_error=np.empty(0)
                r = 1
                l = 2
                while r:
                    str_temp=data_linesplit[i+l].split()
                    if not len(data_linesplit[i+l].split())==0:
                        psip_Er=np.append(psip_Er,float(str_temp[0]))
                        Er=np.append(Er,float(str_temp[1]))
                        Er_error=np.append(Er_error,float(str_temp[2]))
                        l=l+1
                        i=i+1
                    else:
                        r= (1==2)
        i = i+1

        if i == len(data_linesplit):
            keep_going=(1==2)

    return psip_tb,tb,psip_nb,nb,psip_Er,Er

def read_cxrs_file_full(cxrs_filename):

    cxrs_file=open(cxrs_filename,'r')

    data_in=cxrs_file.read()
    data_linesplit=data_in.split('\n')

    keep_going=1
    i=0
    while keep_going:
        test=re.search('psinorm',data_linesplit[i])
        if test:
            quantity=data_linesplit[i].split()[1]
            print "Reading :",quantity
            if quantity=='T_{Boron5+}(keV)':
                psip_tb=np.empty(0)
                tb=np.empty(0)
                tb_error=np.empty(0)
	        r = 1
		l = 2
                while r:
                    str_temp=data_linesplit[i+l].split()
                    if not len(data_linesplit[i+l].split())==0:
                        psip_tb=np.append(psip_tb,float(str_temp[0]))
                        tb=np.append(tb,float(str_temp[1]))
                        tb_error=np.append(tb_error,float(str_temp[2]))
                        l=l+1
                        i=i+1
                    else:
                        r= (1==2)
            if quantity=='n_{Boron5+}(10^18':
                psip_nb=np.empty(0)
                nb=np.empty(0)
                nb_error=np.empty(0)
                r = 1
                l = 2
                while r:
                    str_temp=data_linesplit[i+l].split()
                    if not len(data_linesplit[i+l].split())==0:
                        psip_nb=np.append(psip_nb,float(str_temp[0]))
                        nb=np.append(nb,float(str_temp[1]))
                        nb_error=np.append(nb_error,float(str_temp[2]))
                        l=l+1
                        i=i+1
                    else:
                        r= (1==2)
            if quantity=='vpol_{Boron5+}(km/s)':
                psip_vpolb=np.empty(0)
                vpolb=np.empty(0)
                vpolb_error=np.empty(0)
                r = 1
                l = 2
                while r:
                    str_temp=data_linesplit[i+l].split()
                    if not len(data_linesplit[i+l].split())==0:
                        psip_vpolb=np.append(psip_vpolb,float(str_temp[0]))
                        vpolb=np.append(vpolb,float(str_temp[1]))
                        vpolb_error=np.append(vpolb_error,float(str_temp[2]))
                        l=l+1
                        i=i+1
                    else:
                        r= (1==2)
            if quantity=='vtor_{Boron5+}(km/s)':
                psip_vtorb=np.empty(0)
                vtorb=np.empty(0)
                vtorb_error=np.empty(0)
                r = 1
                l = 2
                while r:
                    str_temp=data_linesplit[i+l].split()
                    if not len(data_linesplit[i+l].split())==0:
                        psip_vtorb=np.append(psip_vtorb,float(str_temp[0]))
                        vtorb=np.append(vtorb,float(str_temp[1]))
                        vtorb_error=np.append(vtorb_error,float(str_temp[2]))
                        l=l+1
                        i=i+1
                    else:
                        r= (1==2)
            if quantity=='Er(kV/m)':
                psip_Er=np.empty(0)
                Er=np.empty(0)
                Er_error=np.empty(0)
                r = 1
                l = 2
                while r:
                    str_temp=data_linesplit[i+l].split()
                    if not len(data_linesplit[i+l].split())==0:
                        psip_Er=np.append(psip_Er,float(str_temp[0]))
                        Er=np.append(Er,float(str_temp[1]))
                        Er_error=np.append(Er_error,float(str_temp[2]))
                        l=l+1
                        i=i+1
                    else:
                        r= (1==2)
        i = i+1

        if i == len(data_linesplit):
            keep_going=(1==2)

    return psip_tb,tb,psip_nb,nb,psip_Er,Er,psip_vpolb,vpolb,psip_vtorb,vtorb
