#! /usr/bin/env python
import sys
import time

def mycall(comand, alias_py='python'):
    alias_py=''
    from subprocess import call
    line=comand.split(" ")
    if alias_py != '' :
        fcomand=[alias_py]
    else:
        fcomand=[]#2.7']
    #fcomand.extend(line)#2.7']
    fcomand.extend(line)
    linp=''
    nx=len(fcomand)
    for i in range(0, nx):
        linp=linp+fcomand[i]+" "
    print linp+" | PROCESSING"
    flog.write(linp+"\n")
    call(fcomand)
    print "DONE"
    
def print_time():
    [year,mon,mday,hour,min,sec,wday,yday,isdst] = time.localtime()
    print "# TIME "+str(sec)+" "+str(min)+" "+str(hour)+" "+str(mday)+" "+str(mon)+" "+str(year)+" "+str(wday)+" "+str(yday)+" "+str(isdst)
    flog.write("# TIME "+str(sec)+" "+str(min)+" "+str(hour)+" "+str(mday)+" "+str(mon)+" "+str(year)+" "+str(wday)+" "+str(yday)+" "+str(isdst)+"\n")
    sec_now=yday*24.*3600.+hour*3600.+min*60.+sec
    return sec_now

def sycall(comand):
    import os
    linp=comand
    print linp+" | PROCESSING"
    flog.write(linp+"\n")
    os.system(comand)
    print "DONE"

sys.argv=filter(None,sys.argv)
if len(sys.argv) < 2:
    print "ana_single.pl NAME [DIR_DATA] [DIR_DATA_OUT] [DIR_PLOTS]";
    sys.exit(0)
DIR_DATA="../disk-a/sanchez/ppak/legacy/DATA/SSP/MaNGA_lin"
DIR_DATA_OUT="../disk-b/manga/data/v1_2_0/MaNGA_lin_ana"
DIR_PLOTS="../home/manga/MaNGA_figs"
NAME=sys.argv[1]
if len(sys.argv) == 3:
    DIR_DATA=sys.argv[2]
if len(sys.argv) == 4:
    DIR_DATA=sys.argv[2]
    DIR_DATA_OUT=sys.argv[3]
if len(sys.argv) == 5:
    DIR_DATA=sys.argv[2]
    DIR_DATA_OUT=sys.argv[3]
    DIR_PLOTS=sys.argv[4]
SPLIT=NAME.split("-",3)
#DIR_DATA_OUT_final0=DIR_DATA_OUT+"/"+SPLIT[1]
DIR_DATA_OUT_final=DIR_DATA_OUT+"/"+NAME#SPLIT[1]+"/"+SPLIT[2]

logfile="ana_single_plot."+NAME+".log"
flog=open(logfile,"w")
init_time=print_time()

DIRS=DIR_DATA_OUT.split("/")
DRT=""
for DR in DIRS:
    DRT=DRT+DR+"/"
    call="mkdir -p "+DRT
    sycall(call) 
#call="mkdir -p "+DIR_DATA_OUT_final0
#sycall(call)
call="mkdir -p "+DIR_DATA_OUT_final
sycall(call)

DIRS=DIR_PLOTS.split("/")
DRT=""
for DR in DIRS:
    DRT=DRT+DR+"/"
    call="mkdir -p "+DRT
    sycall(call) 
call="mkdir -p "+DIR_PLOTS+"/"+SPLIT[1]
sycall(call)
call="mkdir -p "+DIR_PLOTS+"/"+SPLIT[1]+"/"+SPLIT[2]
sycall(call)


call="plot_maps.py "+DIR_DATA_OUT_final+"/map.6530_6630."+NAME+"_flux_00.fits.gz 0 3 1 1 Ha,"+NAME+" 1 "+DIR_PLOTS+"/"+SPLIT[1]+"/"+SPLIT[2]+"/Ha_flux_"+NAME+".pdf"
mycall(call)
call="plot_maps.py "+DIR_DATA_OUT_final+"/"+NAME+".V.fits.gz 0 0.5 1 1 V,"+NAME+" 1 "+DIR_PLOTS+"/"+SPLIT[1]+"/"+SPLIT[2]+"/V_flux_"+NAME+".pdf"
mycall(call)
call="plot_maps.py "+DIR_DATA_OUT_final+"/map.6530_6630."+NAME+"_vel_00.fits.gz 7400 7700 1 1 VEL,"+NAME+" 1 "+DIR_PLOTS+"/"+SPLIT[1]+"/"+SPLIT[2]+"/Ha_vel_"+NAME+".pdf"
mycall(call)

end_time=print_time()
delta_time=end_time-init_time
print "Total Time consumed = "+str(delta_time)

flog.close()
