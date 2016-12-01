#!/usr/bin/python

import sys
import numpy as np
from numpy import float_ as flt
import time
import os.path as pt
import pyfits
from pyfits import getheader as ghead
from pyfits import getdata as gdata
from pyfits import writeto as wdata

def sycallo(comand):
    import os
    out=os.popen(comand, 'r')
    line=out.readline()
    return line.replace("\n","")

def is_mpi():
    import imp
    try:
        imp.find_module('mpi4py')
        found = True
    except ImportError:
        found = False
    return found

def print_time():
    [year,mon,mday,hour,min,sec,wday,yday,isdst] = time.localtime()
    print "# TIME "+str(sec)+" "+str(min)+" "+str(hour)+" "+str(mday)+" "+str(mon)+" "+str(year)+" "+str(wday)+" "+str(yday)+" "+str(isdst)
    flog.write("# TIME "+str(sec)+" "+str(min)+" "+str(hour)+" "+str(mday)+" "+str(mon)+" "+str(year)+" "+str(wday)+" "+str(yday)+" "+str(isdst)+"\n")
    sec_now=yday*24.*3600.+hour*3600.+min*60.+sec
    return sec_now

def mpicall(comand, npr=2, alias_mpi='mpiexec', alias_py='python'):
    alias_py=''
    from subprocess import call
    line=comand.split(" ")
    if alias_py != '':
        fcomand=[alias_mpi+' -np '+str(npr)+' '+alias_py][0].split(" ")
    else:
        fcomand=[alias_mpi+' -np '+str(npr)][0].split(" ")
    #fcomand=['mpiexec2 -np '+str(npr)+' python2.7'][0].split(" ")
    fcomand.extend(line)
    linp=''
    nx=len(fcomand)
    for i in range(3, nx):
        linp=linp+fcomand[i]+" "
    print linp
    flog.write(linp+"\n")
    #print fcomand
    call(fcomand)
    print "DONE"

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
    #print fcomand
    call(fcomand)
    print "DONE"
    
def sycall(comand):
    import os
    linp=comand
    print linp+" | PROCESSING"
    flog.write(linp+"\n")
    os.system(comand)
    print "DONE"
    
def wfits_ext(name, hlist):
    if pt.exists(name) == False:
        hlist.writeto(name)
    else:
        sycall("rm "+name)
        hlist.writeto(name)

def reset_cube(name):
    import pyfits as pyf
    hdu_list=pyf.open(name)
    flux_list= hdu_list['FLUX']
    ivar_list= hdu_list['IVAR']
    wave_list= hdu_list['WAVE']
    mask_list= hdu_list['MASK']
    flux=flux_list.data
    wave=wave_list.data
    nw=len(wave)
    ns=flux.shape[1]
    dl=wave[1]-wave[0]
    h=flux_list.header
    h1_l=pyf.PrimaryHDU(flux)
    h1=h1_l.header
    h1['CRVAL3']=h['CRVAL3']
    h1['CDELT3']=h['CD3_3']
    h1['CRPIX3']=h['CRPIX3']
    h1['CRVAL1']=h['CRVAL1']
    h1['CD1_1']=h['CD1_1']
    h1['CRPIX1']=h['CRPIX1']
    h1['CTYPE1']=h['CTYPE1']
    h1['CRVAL2']=h['CRVAL2']
    h1['CD2_2']=h['CD2_2']
    h1['CRPIX2']=h['CRPIX2']
    h1['CTYPE2']=h['CTYPE2']
    hlist=pyf.HDUList([h1_l,ivar_list,mask_list,mask_list])
    wfits_ext(name.replace('-LINCUBE.fits.gz','.cube.fits'),hlist)
    call="gzip "+name.replace('-LINCUBE.fits.gz','.cube.fits')
    sycall(call)

sys.argv=filter(None,sys.argv)        
if len(sys.argv) < 1:
    print "USE: ana_MaNGA.py NAME"
    sys.exit(0)
NAME=sys.argv[1]
logfile="ana_MaNGA."+NAME.replace('-LINCUBE.fits.gz','')+".log"
flog=open(logfile,"w")
if pt.exists(NAME.replace('-LINCUBE.fits.gz','.cube.fits.gz')) == False:
    reset_cube(NAME)
#print NAME
#sys.exit()
NAME=NAME.replace("-LINCUBE.fits.gz","")

n_proc=1
vel_light=flt(299792.458)
plot=2 
plot_rss=0

VER="v1_3_2"
VER="MPL-5"
DIR_DATA="./"
DIR_DATA_OUT="./"
DIR_PLOTS="./"
max_size=10
frac=0.9
#template="home/sanchez/ppak/legacy/gsd61_156.fits"
#temp_2="disk-b/sanchez/ppak/legacy/miles_2.fits"
#temp_3="home/sanchez/ppak/legacy/gsd61_12.fits"
#temp_4="disk-b/sanchez/ppak/legacy/miles_2_gas.fits"
#temp_5="disk-b/sanchez/ppak/legacy/templates/ssp_lib.4.fits"
#temp_6="disk-b/sanchez/ppak/legacy/templates/ssp_lib.24.fits"
DIR_CONF="legacy"
template="disk-b/sanchez/ppak/legacy/gsd01_156.fits"
temp_2="disk-b/sanchez/ppak/legacy/miles_2.fits"
temp_3="disk-b/sanchez/ppak/legacy/templates/miuscat/spectra/miuscat_12.fits"
temp_4="disk-b/sanchez/ppak/legacy/miles_2_gas.fits"
temp_5="disk-b/sanchez/ppak/legacy/templates/miuscat/spectra/miuscat_12.fits"
temp_6="disk-b/sanchez/ppak/legacy/templates/miuscat/spectra/miuscat_12.fits"
pack_type="v1.5"
config="auto_ssp_MaNGA_strong.config"
config_no="auto_ssp_MaNGA_no_lines.config"
config_SII="auto_ssp_MaNGA_strong_SII.config"

root_dir=sycallo("echo $FIT3DP_PATH")+"/.."
DIR_CONF=root_dir+"/"+DIR_CONF
#DIR_DATA=root_dir+"/"+DIR_DATA
#DIR_DATA_OUT=root_dir+"/"+DIR_DATA_OUT
#DIR_PLOTS=root_dir+"/"+DIR_PLOTS
template=root_dir+"/"+template
temp_2=root_dir+"/"+temp_2
temp_3=root_dir+"/"+temp_3
temp_4=root_dir+"/"+temp_4
temp_5=root_dir+"/"+temp_5
temp_6=root_dir+"/"+temp_6
nsa_z=0;
if len(sys.argv) == 3:
    nsa_z=sys.argv[2]
SPLIT=NAME.split("-",3)
#DIR_DATA_OUT_final0=DIR_DATA_OUT+"/"+SPLIT[1]
DIR_DATA_OUT_final=DIR_DATA_OUT+"/"+NAME#+"/"+SPLIT[2]
FILE=DIR_DATA+NAME+".cube.fits.gz"
init_time=print_time()
call="re_config.py "+config+" "+config_no+" "+config_SII+" "+DIR_CONF
mycall(call)
#sys.exit()
config=config+"_n"
config_no=config_no+"_n"
config_SII=config_SII+"_n"
call="get_slice.py "+FILE+" img_"+NAME+" "+DIR_CONF+"/slice_V.conf" 
mycall(call)

V_img="img_"+NAME+"_V_4500_5500.fits"
[img,hd_mask]=gdata(V_img,0,header=True)                                             
[nx,ny]=img.shape    
img_mask=np.zeros([nx,ny])
V_mask="../masks_v1.4/"+NAME+".mask.fits"
nx1=nx 
ny1=ny
if pt.exists(V_mask) == True:
    [img_mask,hd_mask]=gdata(V_mask,0,header=True)
    [nx1,ny1]=img_mask.shape   
print "DIM = "+str(nx)+","+str(ny)+","+str(nx1)+","+str(ny1)          
img_mask_new=np.ones([nx1,ny1])
for ii in range(0, nx1):                                                     
    for jj in range(0, ny1):
        val_mask=img_mask[ii,jj]                                                 
        if val_mask == 2:
            img_mask_new[ii,jj]=0
        else:
            img_mask_new[ii,jj]=1            
#
# We look for stars in the field
#
ic=int(nx*0.5)
jc=int(ny*0.5)
val_c=img[ic,jc]
dist_max=0.15*np.sqrt(nx**2+ny**2)
for ii in range(0, nx):
    for jj in range(0, ny):     
        val=img[ii,jj]
        dist=np.sqrt((ii-ic)**2+(jj-jc)**2)
        if ((val > val_c) and (dist > dist_max)):
            # We mask the star
            nx0=ii-3
            nx1=ii+3
            ny0=jj-3
            ny1=jj+3
            if nx0 < 0:
                nx0=0
            if nx1 > (nx-1):
                nx1=nx-1
            if ny0 < 0:
                ny0=0
            if ny1 > (ny-1):
                ny1=ny-1
            for ii2 in range(nx0, nx1):                                                  
                for jj2 in range(ny0, ny1):                                    
                    img_mask_new[ii2,jj2]=1
                    
V_mask=NAME+".mask.fits"
if pt.exists(V_mask) == False:
    wdata(V_mask,img_mask_new,hd_mask)
else:
    sycall("rm "+V_mask)
    wdata(V_mask,img_mask_new,hd_mask)

print "PASO"                        
val_max=-1e12
XC=0                                                                
YC=0                        
SUM=0                                                
nx0=int(nx*0.33)
nx1=int(nx*2*0.33)
ny0=int(nx*0.33)
ny1=int(nx*2*0.33)
for ii in range(nx0, nx1):                                                  
    for jj in range(ny0, ny1):                                                 
        val=img[ii,jj]                                                 
        val_mask=img_mask[ii,jj]                                                 
        if val_mask != 2:
            XC=XC+ii*(val)**5                                                            
            YC=YC+jj*(val)**5                                                           
            SUM=SUM+(val)**5   
XC=XC/SUM
YC=YC/SUM
call="cp "+V_img+" "+NAME+".V.fits"
sycall(call)

call="clean_nan.py "+NAME+".V.fits -1"
mycall(call)

call="radial_sum_cube_e.py "+FILE+" 2.5 "+str(int((np.round(XC))))+" "+str(int((np.round(YC))))+" rad."+NAME+".rss.fits 2"
mycall(call)

call="img2spec_e.py rad."+NAME+".rss.fits 0 "+NAME+".spec_5.txt"
mycall(call)

call="img2spec_e.py rad."+NAME+".rss.fits 5 "+NAME+".spec_30.txt"
mycall(call)

call="read_img_header.py "+FILE+" MED_VEL SYS_VEL_NOW_"+NAME
mycall(call)

sys_vel_now="SYS_VEL_NOW_"+NAME
f=open(sys_vel_now,'r')
for line in f:
    data=line.split(" ")
f.close()
key=data[1]
vel=data[2]
if vel== "none":
    vel=0
call="vel_eline_spec.py "+NAME+".spec_5.txt 2 0.7 junk 3728 OII_"+NAME+".pdf/1 3700 4000 vel_eline_spec."+NAME+".txt"
mycall(call)

f_eline="vel_eline_spec."+NAME+".txt"
f=open(f_eline,'r')
for line in f:
    dat_fe=line.split(",")
f.close()
#
# Define the redshift by OII, if it is clearly detected!
#
if dat_fe[1] == '1' and vel== 0:
    vel=flt(dat_fe[2])
    red=vel/vel_light;
    min_red=red-300/vel_light;
    max_red=red+300/vel_light;
    d_red=30/vel_light
else:
    red=0.02
    min_red=0.0045
    max_red=0.032
    d_red=30/vel_light
red=0.02
min_red=0.0001
max_red=0.15
d_red=300/vel_light
call="auto_ssp_elines_rnd.py  "+NAME+".spec_5.txt "+temp_5+","+temp_4+" auto_ssp_Z."+NAME+".cen.out "+DIR_CONF+"/mask_elines.txt "+DIR_CONF+"/"+config_no+" "+str(plot)+" -3 50 3800,3850 7000,4700 none  "+str(red)+" "+str(d_red)+" "+str(min_red)+" "+str(max_red)+" 1.6 0.5 1.2 5.5 0.0 0.15 0.0 1.6"
mycall(call)
#sys.exit()
vel=0
disp_MAX=5.5
if vel==0:
    outfile="auto_ssp_Z."+NAME+".cen.out"
    f=open(outfile, 'r')
    for line in f:
        RED_DATA=line.split(",")
        if not '#' in RED_DATA[0]:
            red=flt(RED_DATA[7])
            min_red=red-300/vel_light
            max_red=red+300/vel_light
            d_red=30/vel_light
            vel=red*vel_light
            disp_MAX=1.5*flt(RED_DATA[9])
            print "RED="+str(red)+ "VEL="+str(vel)
    f.close()

call="redshift_config.py "+DIR_CONF+"/"+config+" "+str(red)+" auto."+NAME+".config"
mycall(call)

w_cut1_1=3800*(1+red)
w_cut1_2=7000*(1+red)
w_cut2_1=3850*(1+red)
w_cut2_2=4700*(1+red)

#call="auto_ssp_elines_rnd.py  "+NAME+".spec_5.txt "+temp_5+","+temp_2+" auto_ssp_no_mask."+NAME+".cen.out none auto."+NAME+".config "+str(plot)+" -3 50 3700,3850 10000,4700 "+DIR_CONF+"/emission_lines.txt  "+str(red)+" 0 "+str(min_red)+" "+str(max_red)+" 3.2 0.3 1.2 "+str(disp_MAX)+" 0.3 0.15 0.0 1.6"
call="auto_ssp_elines_rnd_sigma_inst.py  "+NAME+".spec_5.txt "+temp_5+","+temp_2+",0.7201 auto_ssp_no_mask."+NAME+".cen.out none auto."+NAME+".config "+str(plot)+" -3 50 "+str(w_cut1_1)+","+str(w_cut2_1)+" "+str(w_cut1_2)+","+str(w_cut2_2)+" "+DIR_CONF+"/emission_lines.txt  "+str(red)+" 0 "+str(min_red)+" "+str(max_red)+" 30 20 1 400 0.0 0.15 0.0 1.6"
mycall(call)

#call="auto_ssp_elines_rnd.py  "+NAME+".spec_5.txt "+temp_5+","+temp_2+" auto_ssp."+NAME+".cen.out "+DIR_CONF+"/mask_elines.txt auto."+NAME+".config "+str(plot)+" -3 50 3700,3850 10000,4700 "+DIR_CONF+"/emission_lines.txt  "+str(red)+" 0 "+str(min_red)+" "+str(max_red)+" 3.2 0.3 1.2 "+str(disp_MAX)+" 0.3 0.15 0.0 1.6"
call="auto_ssp_elines_rnd_sigma_inst.py  "+NAME+".spec_5.txt "+temp_5+","+temp_2+",0.7201 auto_ssp."+NAME+".cen.out "+DIR_CONF+"/mask_elines.txt auto."+NAME+".config "+str(plot)+" -3 50 "+str(w_cut1_1)+","+str(w_cut2_1)+" "+str(w_cut1_2)+","+str(w_cut2_2)+" "+str(DIR_CONF)+"/emission_lines.txt  "+str(red)+" 0 "+str(min_red)+" "+str(max_red)+" 30 20 1 400 0.0 0.15 0.0 1.6"
mycall(call)

#call="auto_ssp_elines_rnd.py  "+NAME+".spec_30.txt "+temp_5+","+temp_2+" auto_ssp."+NAME+".int.out "+DIR_CONF+"/mask_elines.txt auto."+NAME+".config "+str(plot)+" -3 150 3700,3850 10000,4700 "+DIR_CONF+"/emission_lines.txt  "+str(red)+" "+str(d_red)+" "+str(min_red)+" "+str(max_red)+" 3.2 0.3 1.2 "+str(disp_MAX)+" 0.3 0.15 0.0 1.6"
call="auto_ssp_elines_rnd_sigma_inst.py  "+NAME+".spec_30.txt "+temp_5+","+temp_2+",0.7201 auto_ssp."+NAME+".int.out "+DIR_CONF+"/mask_elines.txt auto."+NAME+".config "+str(plot)+" -3 150 "+str(w_cut1_1)+","+str(w_cut2_1)+" "+str(w_cut1_2)+","+str(w_cut2_2)+" "+DIR_CONF+"/emission_lines.txt  "+str(red)+" "+str(d_red)+" "+str(min_red)+" "+str(max_red)+" 30 20 1 400 0.0 0.15 0.0 1.6"
mycall(call)
#sys.exit()
call="med2df.py "+NAME+".V.fits m"+NAME+".V.fits 2 2"
mycall(call)

call="imarith.py m"+NAME+".V.fits '*' "+V_mask+" m"+NAME+".V.fits"
mycall(call)

call="imarith.py m"+NAME+".V.fits '*' "+V_mask+" "+NAME+".V.fits"
mycall(call)

call="get_SN_cube.py "+FILE+" map_SN."+NAME+".fits signal."+NAME+".fits noise."+NAME+".fits"
mycall(call)

call="cont_seg_all_SN.py signal."+NAME+".fits,noise."+NAME+".fits 0.001 50,1 0.85 0.001 cont_seg."+NAME+".fits DMASK."+NAME+".fits"
#call="cont_seg_all.py "+NAME+".V.fits 0.001 4 0.85 0.001 cont_seg."+NAME+".fits DMASK."+NAME+".fits"
mycall(call)

#call="create_mask_map.py m"+NAME+".V.fits 0.005 mask."+NAME+".V.fits"
call="create_mask_map.py map_SN."+NAME+".fits 1 mask."+NAME+".V.fits"
mycall(call)

call="spec_extract_cube_mean.py "+FILE+" cont_seg."+NAME+".fits CS."+NAME+".RSS.fits"
mycall(call)

call="spec_extract_cube_error.py "+FILE+" cont_seg."+NAME+".fits e_CS."+NAME+".RSS.fits"
mycall(call)

call="get_SN_rss.py CS."+NAME+".RSS.fits SN_"+NAME+".CS.rss.csv"
mycall(call)

call="csv_to_map.seg.py SN_"+NAME+".CS.rss.csv 1 cont_seg."+NAME+".fits SN_"+NAME+".CS.fits"
mycall(call)

call="rss_seg2cube.py CS."+NAME+".RSS.fits cont_seg."+NAME+".fits SEG.cube."+NAME+".fits"
mycall(call)

call="get_slice.py SEG.cube."+NAME+".fits SEG_img_"+NAME+" "+DIR_CONF+"/slice_V.conf" 
mycall(call)

call="clean_nan.py  SEG_img_"+NAME+"_V_4500_5500.fits -1"
mycall(call)

call="imarith.py "+NAME+".V.fits '/' SEG_img_"+NAME+"_V_4500_5500.fits scale.seg."+NAME+".fits"
mycall(call)

print "Fitting RSS file"

#Check if MPI is installed
f_mpi=is_mpi()
if f_mpi:
    if n_proc > 1:
        call_mpi="_mpi"
        print "MPI4PY MODULE FOUND. FITTTING RSS FILE IN PARALLEL MODE"
    else:
        call_mpi=""
else:
    if n_proc > 1:
        print "WARNING: THERE IS NO MPI4PY MODULE INSTALLED, FITTTING RSS FILE IN NON-PARALLEL MODE"
    call_mpi=""
call="auto_ssp_elines_rnd_rss_sigma_inst_disp_min"+call_mpi+".py  CS."+NAME+".RSS.fits,e_CS."+NAME+".RSS.fits "+temp_5+","+temp_2+",0.7201  auto_ssp.CS_few."+NAME+".rss.out "+DIR_CONF+"/mask_elines.txt auto."+NAME+".config "+str(plot_rss)+" -2 5 "+str(w_cut1_1)+","+str(w_cut2_1)+" "+str(w_cut1_2)+","+str(w_cut2_2)+" "+DIR_CONF+"/emission_lines.txt  "+str(red)+" "+str(d_red)+" "+str(min_red)+" "+str(max_red)+"  30 20 1 400 0.0 0.15 0.0 1.6"
#call="auto_ssp_elines_rnd_rss_sigma_inst_disp_min"+call_mpi+".py  CS."+NAME+".RSS.fits,e_CS."+NAME+".RSS.fits "+temp_5+","+temp_2+" auto_ssp.CS_few."+NAME+".rss.out "+DIR_CONF+"/mask_elines.txt "+DIR_CONF+"/"+config_no+" "+str(plot_rss)+" -2 5 3700,3850 10000,4700 "+DIR_CONF+"/emission_lines.txt  "+str(red)+" "+str(d_red)+" "+str(min_red)+" "+str(max_red)+" 3.2 0.3 1.2 "+str(disp_MAX)+" 0.3 0.15 0.0 1.6"
if f_mpi:
    if n_proc > 1:
        mpicall(call,n_proc)
    else:
        mycall(call)
else:
    mycall(call)
    
#sys.exit()
    
call="plot_output_auto_ssp_elines_several_Av_log_rss_all.py output.auto_ssp.CS_few."+NAME+".rss.out.fits plot_fit.CS_few."+NAME+".pdf"
mycall(call)

call="FIT3D_output_rss_seg2cube.py output.auto_ssp.CS_few."+NAME+".rss.out.fits 1 cont_seg."+NAME+".fits SSP_mod_tmp."+NAME+".cube.fits"
mycall(call)

call="imarith.py SSP_mod_tmp."+NAME+".cube.fits '*' scale.seg."+NAME+".fits SSP_mod."+NAME+".cube.fits"
mycall(call)

call="imarith.py "+FILE+" '-' SSP_mod."+NAME+".cube.fits TMP."+NAME+".cube.fits"
mycall(call)

call="rm -r smooth."+NAME+".cube.fits"
sycall(call)

call="smooth_spec_clip_cube.py TMP."+NAME+".cube.fits smooth."+NAME+".cube.fits 75 1.5 10 1860"
mycall(call)

call="imarith.py TMP."+NAME+".cube.fits - smooth."+NAME+".cube.fits GAS."+NAME+".cube.fits"
mycall(call)

print " GAS."+NAME+".cube.fits created"

call="spec_extract_cube_mean.py GAS."+NAME+".cube.fits cont_seg."+NAME+".fits GAS.CS."+NAME+".RSS.fits"
mycall(call)

disp_MAX=4.7
call="auto_ssp_elines_rnd_rss_sigma_inst_guided"+call_mpi+".py  output.auto_ssp.CS_few."+NAME+".rss.out.fits,e_CS."+NAME+".RSS.fits "+template+","+temp_2+",auto_ssp.CS_few."+NAME+".rss.out,1.601 auto_ssp.CS."+NAME+".rss.out "+DIR_CONF+"/mask_elines.txt "+DIR_CONF+"/"+config_no+" "+str(plot_rss)+" -2 5 "+str(w_cut1_1)+","+str(w_cut2_1)+" "+str(w_cut1_2)+","+str(w_cut2_2)+" "+DIR_CONF+"/emission_lines.txt  "+str(red)+" "+str(d_red)+" "+str(min_red)+" "+str(max_red)+" 30 20 1 400 0.0 0.15 0.0 1.6"
#call="auto_ssp_elines_rnd_rss_disp_MaNGA"+call_mpi+".py  CS."+NAME+".RSS.fits,e_CS."+NAME+".RSS.fits "+template+","+temp_2+" auto_ssp.CS."+NAME+".rss.out "+DIR_CONF+"/mask_elines.txt auto."+NAME+".config "+str(plot_rss)+" -2 5 3800,3850 7000,4700 "+DIR_CONF+"/emission_lines.txt  "+str(red)+" "+str(d_red)+" "+str(min_red)+" "+str(max_red)+" 2.5 0.15 0.6 "+str(disp_MAX)+" 0.3 0.15 0.0 1.6"
if f_mpi:
    if n_proc > 1:
        mpicall(call,n_proc)
    else:
        mycall(call)
else:
    mycall(call)
call="cp elines_auto_ssp.CS_few."+NAME+".rss.out elines_auto_ssp.CS."+NAME+".rss.out"
sycall(call)

call="plot_output_auto_ssp_elines_several_Av_log_rss_all.py output.auto_ssp.CS."+NAME+".rss.out.fits plot_fit.CS."+NAME+".pdf"
mycall(call)
#
# We analyze the indeces.
#
call="get_index_output_auto_ssp_elines_rnd_rss_outFIT3D.py output.auto_ssp.CS."+NAME+".rss.out.fits 5 auto_ssp.CS."+NAME+".rss.out /null indices.CS."+NAME+".rss.out "
mycall(call)
#sys.exit()
call="index_seg_cube.py indices.CS."+NAME+".rss.out cont_seg."+NAME+".fits indices.CS."+NAME+".cube.fits"
mycall(call)
#sys.exit()
call="map_auto_ssp_rnd_seg.py elines_auto_ssp.CS."+NAME+".rss.out cont_seg."+NAME+".fits map.CS."+NAME+" 0.9 "
mycall(call)
#sys.exit()
WREF=6562.68
WMIN=WREF*(1+min_red)-50
WMAX=WREF*(1+max_red)+50
call="vel_eline_cube.py GAS."+NAME+".cube.fits 1 0.95 ve."+NAME+" "+str(WREF)+" /null "+str(WMIN)+" "+str(WMAX)
mycall(call)

min_vel=red*vel_light-500
max_vel=red*vel_light+500
call="clean_Ha_map.py   ve."+NAME+".vel_map.fits ve."+NAME+".mask_map.fits "+str(max_vel)+" "+str(min_vel)
mycall(call)
#
# Analyze cube
#

elines_in_file="elines_auto_ssp.CS."+NAME+".rss.out"
f=open(elines_in_file,"r")
line=f.readline()
line=f.readline()
line=f.readline()
data_elines=line.split(" ")
data_elines=filter(None,data_elines)
f.close
disp_max_elines=flt(data_elines[5])*1.05

FILE_CONF_NOW=DIR_CONF+"/"+config_SII
f=open(FILE_CONF_NOW, 'r')
conf=f.readline()
conf=f.readline()
conf=f.readline()
nr=f.readline()
nr=nr.split(" ")
nr=int(filter(None,nr)[0])
print "NR="+str(nr)
for j in range(0, nr):
    conf=f.readline()
    d_conf=conf.split(" ")
    d_conf=filter(None,d_conf)
    w1=flt(d_conf[0])
    w2=flt(d_conf[1])
    w1_z=int(w1*(1+red))
    w2_z=int(w2*(1+red))
    in_conf=d_conf[3]
    f2=open(in_conf, 'r')
    tmp_conf="tmp."+NAME+".conf"
    f3=open(tmp_conf, 'w')
    ni=0
    for line in f2:
        ins=line
        if ni == 5:
            vel1=vel-300
            vel2=vel+300
            f3.write(str(vel)+"     1       "+str(vel1)+"     "+str(vel2)+"   -1\n")
        else:
            if ni == 4:
                f3.write("1.2     1       0.75     "+str(disp_max_elines)+"   -1\n")
            else:
                f3.write(ins)
        ni=ni+1
    f2.close()
    f3.close()
    print w1,w2
    if j == 0:
    # We found not difference in the Kinematics fitted and guessed with an hyperbolic curve
    # Thus, we fixed the velocity map: 
        call="kin_cube_elines_rnd.py  GAS."+NAME+".cube.fits  "+tmp_conf+" none "+str(w1_z)+" "+str(w2_z)+" KIN.GAS."+str(int(w1))+"_"+str(int(w2))+"."+NAME+".out 30 3 0 0.15 map."+str(int(w1))+"_"+str(int(w2))+"."+NAME+" 0 ve."+NAME+".vel_map.fits,ve."+NAME+".mask_map.fits,0 "#>/dev/null "
    else:
        fix_disp=0
        if w1 == 3700:
            fix_disp=1
        call="kin_cube_elines_rnd.py  GAS."+NAME+".cube.fits  "+tmp_conf+" none "+str(w1_z)+" "+str(w2_z)+" KIN.GAS."+str(int(w1))+"_"+str(int(w2))+"."+NAME+".out 30 3 0 0.15 map."+str(int(w1))+"_"+str(int(w2))+"."+NAME+" 0 map.6530_6630."+NAME+"_vel_00.fits,ve."+NAME+".mask_map.fits,0 map.6530_6630."+NAME+"_disp_00.fits,"+str(fix_disp)+" "#>/dev/null"
    mycall(call)
    if j == 0:
        call="rm -f vel.fits"
        sycall(call)
        call="rm -f disp.fits"
        sycall(call)
        call="med2df.py map.6530_6630."+NAME+"_vel_00.fits vel."+NAME+".fits 3 3"
        mycall(call);
        call="med2df.py map.6530_6630."+NAME+"_disp_00.fits disp."+NAME+".fits 3 3"
        mycall(call)
        call="imarith.py disp."+NAME+".fits '/' 2.345 disp."+NAME+".fits"
        mycall(call)
    print str(w1)+"-"+str(w2)+" region analyzed"
f.close
#
# HII regions
#

call="gzip -f *"+NAME+"*.fits"
sycall(call)

call="get_CS_slice.py "+NAME
mycall(call)

call="ana_single_elines.py "+NAME+" "+DIR_DATA+" "+DIR_CONF
mycall(call)

call="sum_mass_age.py "+NAME+" MASS."+NAME+".csv"
mycall(call)
#
call="pack_NAME.py "+NAME+" "+DIR_CONF+" "+pack_type
mycall(call)
#sys.exit()
end_time=print_time()
delta_time=end_time-init_time
print "Total Time consumed = "+str(delta_time)
#
# COPY RESULTS!
#
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
call="mv *"+NAME+"* "+DIR_DATA_OUT_final+"/"
sycall(call)
call="mv *.pdf "+DIR_DATA_OUT_final+"/"
sycall(call)
call="mv *.fits "+DIR_DATA_OUT_final+"/"
sycall(call)

call="ana_single_MaNGA_plot_2.py "+NAME+" "+DIR_DATA+" "+DIR_DATA_OUT+" "+DIR_PLOTS
mycall(call)
flog.close
call="mv *.log "+DIR_DATA_OUT_final+"/"
sycall(call)
