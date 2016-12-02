#! /usr/bin/env python
import matplotlib
import sys
import numpy as np
from pyfits import getheader as ghead
from pyfits import getdata as gdata
from pyfits import writeto as wfits
import pyfits as pyf
from scipy.interpolate.interpolate import interp1d

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

sys.argv=filter(None,sys.argv) 

root_dir=sycallo("echo $FIT3DP_PATH")+"/../"
name=sys.argv[1].replace('.fits','')
namef=name+".fits"
hdu_list = pyf.open(namef)
table_hdu = hdu_list[1]
table_hdu1 = hdu_list[2]
table_data = table_hdu.data
table_data1 = table_hdu1.data
flux=table_data.field('flux')
flux_e=table_data.field('ivar')
logl=table_data.field('loglam')
#red=table_data1.field('Z')[0]
wave_min=3800.0
wave_max=9100.0
dl=1.25
wave=np.arange(wave_min,wave_max,dl)
fluxl = interp1d(10**logl, flux/10.0,kind='linear',bounds_error=False)(wave)
fluxl_e = interp1d(10**logl, flux_e,kind='linear',bounds_error=False)(wave)
f=open(name+".spec.txt","w")
for i in range(0, len(wave)):
    if np.isnan(fluxl[i]) == False:
        f.write(str(i+1)+"  "+str(wave[i])+"  "+str(fluxl[i])+"  "+str(fluxl_e[i])+"\n")
    else:
        f.write(str(i+1)+"  "+str(wave[i])+"  "+str(0)+"  "+str(1)+"\n")
f.close()

disp_i=0.21
dirlib=  root_dir+"disk-b/sanchez/ppak/legacy/gsd01_156.fits"
dirlib2= root_dir+"disk-b/sanchez/ppak/legacy/gsd01_3.fits"
dirlib3= root_dir+"disk-b/sanchez/ppak/legacy/templates/miuscat/spectra/miuscat_12.fits"
dirlib4= root_dir+"disk-b/sanchez/ppak/legacy/miles_2_gas.fits"
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

#ESTIMATION OF THE REDSHIFT OF THE GALAXY
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
    
#Redshifting the emision lines to be fit
call="redshift_config.py "+linc2+" "+str(red)+" auto."+name+".config" 
sycall(call)

w_cut1_1=4000*(1+red)
w_cut1_2=7000*(1+red)
w_cut2_1=4050*(1+red)
w_cut2_2=5800*(1+red)
#Fiting the spectrum
call="auto_ssp_elines_rnd_sigma_inst.py "+name+".spec.txt "+dirlib+","+dirlib2+","+str(disp_i)+" auto_ssp."+name+".out "+maskf+" auto."+name+".config 2 -1 25 "+str(w_cut1_1)+","+str(w_cut2_1)+" "+str(w_cut1_2)+","+str(w_cut2_2)+" "+emlines+" "+str(red)+" 0.0001 "+str(red-0.001)+" "+str(red+0.002)+" 50 5 30 300 0.5 0.1 0.0 1.6"
sycall(call)
#Saving the Outputs
call="mkdir -p "+name
sycall(call)
call="mv *.pdf "+name+"/"
sycall(call)
call="mv *.out "+name+"/"
sycall(call)
call="mv *.txt "+name+"/"
sycall(call)
call="mv output.auto_ssp."+name+".out.fits "+name+"/"
sycall(call)