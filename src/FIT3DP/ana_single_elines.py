#!/usr/bin/python
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

sys.argv=filter(None,sys.argv)
if len(sys.argv) < 2:
    print "ana_single_elines.pl NAME [DIR_DATA] [DIR_CONF]"
    sys.exit(0)
NAME=sys.argv[1]
DIR_DATA="../disk-d/sanchez/ppak/legacy/DATA/SSP/MaNGA_lin"
DIR_CONF="../legacy"
if len(sys.argv) == 3:
    DIR_DATA=sys.argv[2]
if len(sys.argv) == 4:
    DIR_DATA=sys.argv[2]
    DIR_CONF=sys.argv[3]
logfile="ana_single_elines."+NAME+".log"
flog=open(logfile,"w")
init_time=print_time()
call="flux_elines_cube_EW.py GAS."+NAME+".cube.fits.gz,"+DIR_DATA+"/"+NAME+".cube.fits.gz[1],10,SSP_mod_tmp."+NAME+".cube.fits.gz "+DIR_CONF+"/emission_lines.MaNGA flux_elines."+NAME+".cube.fits.gz map.6530_6630."+NAME+"_vel_00.fits.gz map.6530_6630."+NAME+"_disp_00.fits.gz"
mycall(call)   
end_time=print_time()
delta_time=end_time-init_time
print "Total Time consumed = "+str(delta_time)
flog.write("Total Time consumed = "+str(delta_time)+" \n")
flog.close()