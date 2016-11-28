import matplotlib
import sys
import numpy as np
from pyfits import getheader as ghead
from pyfits import getdata as gdata
from pyfits import writeto as wfit
import pyfits as pyf
from scipy.interpolate.interpolate import interp1d
import os.path as pt

def sycall(comand):
    print comand
    import os
    linp=comand
    os.system(comand)

def sycallo(comand):
    import os
    out=os.popen(comand, 'r')
    line=out.readline()
    return line.replace("\n","")

def wfits(name, data, hdr):
    if pt.exists(name) == False:
        wfit(name,data,hdr)
    else:
        name1=name.replace("\ "," ")
        name1=name1.replace(" ","\ ")
        sycall("rm "+name1)
        wfit(name,data,hdr)

sys.argv=filter(None,sys.argv) 
root_dir=sycallo("echo $FIT3DP_PATH")+"/../"
name=sys.argv[1]

disp_i=0.7201
dirlib=  root_dir+"disk-b/sanchez/ppak/legacy/gsd01_156.fits"
dirlib2= root_dir+"disk-b/sanchez/ppak/legacy/gsd01_3.fits"
dirlib2= root_dir+"disk-b/sanchez/ppak/legacy/miles_2.fits"
dirlib3= root_dir+"disk-b/sanchez/ppak/legacy/templates/miuscat/spectra/miuscat_12.fits"
dirlib4= root_dir+"disk-b/sanchez/ppak/legacy/miles_2_gas.fits"
#temp_4=  root_dir+"disk-b/sanchez/ppak/legacy/miles_2_gas.fits"
maskf=   root_dir+"legacy/mask_elines.txt"
linc=    root_dir+"legacy/auto_ssp_V500_several_Hb.config"
emlines= root_dir+"legacy/emission_lines.txt"
config="auto_ssp_example_strong.config"
config_no="auto_ssp_MaNGA_no_lines.config"
config_SII="auto_ssp_MaNGA_strong_SII.config"
call="re_config_single.py "+config+" ./"
sycall(call)
call="re_config_single.py "+config_no+" "+root_dir+"legacy/"
sycall(call)
call="re_config_single.py "+config_SII+" "+root_dir+"legacy/"
sycall(call)
config=config+"_n"
config_no=config_no+"_n"
config_SII=config_SII+"_n"
linc=    root_dir+"legacy/"+config_no
linc2=   config

call="rss2spec_e.py CS."+name+".RSS.fits 0 "+name+".spec.txt"
sycall(call)

vel_light=np.float(299792.458)

red=0.01
min_red=0.0001
max_red=0.15
d_red=300/vel_light


call="auto_ssp_elines_rnd.py "+name+".spec.txt "+dirlib3+","+dirlib4+" auto_ssp_Z."+name+".out "+maskf+" "+linc+" 2 -3 50 3800,3850 7000,5700 none  "+str(red)+" "+str(d_red)+" "+str(min_red)+" "+str(max_red)+" 1.6 0.5 1.2 5.5 0.0 0.15 0.0 1.6"
sycall(call)

vel=0
disp_MAX=5.5
if vel==0:
    outfile="auto_ssp_Z."+name+".out"
    f=open(outfile, 'r')
    for line in f:
        RED_DATA=line.split(",")
        if not '#' in RED_DATA[0]:
            red=np.float(RED_DATA[7])
            min_red=red-300/vel_light
            max_red=red+300/vel_light
            d_red=30/vel_light
            vel=red*vel_light
            disp_MAX=1.5*np.float(RED_DATA[9])
            print "RED="+str(red)+ "VEL="+str(vel)
    f.close()

call="redshift_config.py "+linc2+" "+str(red)+" auto."+name+".config"
sycall(call)


w_cut1_1=4000*(1+red)
w_cut1_2=7000*(1+red)
w_cut2_1=4000*(1+red)
w_cut2_2=6800*(1+red)


call="auto_ssp_elines_rnd_rss_sigma_inst_disp_min.py  CS."+name+".RSS.fits,e_CS."+name+".RSS.fits "+dirlib+","+dirlib2+","+str(disp_i)+" auto_ssp.CS."+name+".rss.out "+maskf+" auto."+name+".config 0 -2 5 "+str(w_cut1_1)+","+str(w_cut2_1)+" "+str(w_cut1_2)+","+str(w_cut2_2)+" "+emlines+"  "+str(red)+" "+str(d_red)+" "+str(min_red)+" "+str(max_red)+"  30 20 1 400 0.0 0.15 0.0 1.6"
#print call
sycall(call)
call="plot_output_auto_ssp_elines_several_Av_log_rss_all.py output.auto_ssp.CS."+name+".rss.out.fits plot_fit.CS."+name+".pdf"
sycall(call)
sys.exit()
call="mkdir -p "+name
sycall(call)
call="mv *"+name+"*.pdf "+name+"/"
sycall(call)
call="mv *"+name+"*.out "+name+"/"
sycall(call)
call="mv *"+name+"*.fits.gz "+name+"/"
sycall(call)
call="mv *"+name+"*.RSS.fits "+name+"/"
sycall(call)
call="mv *"+name+"*.txt "+name+"/"
sycall(call)
call="mv *"+name+"*.out.fits "+name+"/"
sycall(call)
call="mv *"+name+"*.config "+name+"/"
sycall(call)