#! /usr/bin/env python
import sys
import numpy as np
from numpy import float_
from numpy import absolute as abs
from numpy import random as ran
import time
import my_auto_ssp_elines_rnd as ssp
import pyfits as pyf
from pyfits import getheader as ghead
from pyfits import getdata as gdata
from pyfits import writeto as wfits
from scipy.interpolate.interpolate import interp1d
import my as my
import os.path as pt
import matplotlib

def sycall(comand):
    from subprocess import call
    line=comand.split(" ")
    fcomand=[]
    fcomand.extend(line)
    linp=''
    nx=len(fcomand)
    for i in range(0, nx):
        linp=linp+fcomand[i]+" "
    print linp
    call(fcomand)

sys.argv=filter(None,sys.argv)
ran.seed(None)
vel_light=299792.458
red_elines=0.0
sec_ini=ssp.print_time()
time1=ssp.get_time()
if len(sys.argv) < 7:
    print "USE: auto_ssp.py SPEC1.txt SSP_SFH.fits,SSP_KIN.fits OUTFILE MASK_LIST CONFIG_FILE PLOT [min max] [wmin wmax] [redshift_elines_to_mask] [input_redshift delta_redshift min_redshift max_redshift] [input_sigma delta_sigma min_sigma max_sigma] [input_Av delta_Av min_Av max_Av]"
    print "CONFIG_FILE:"
    print "redshift delta_redshift min_redshift max_redshift"
    print "sigma delta_sigma min_sigma max_sigma"
    print "Av delta_Av min_Av max_Av [Same range for all]"
    print "N_SYSTEMS"
    print "(1) START_W END_W MASK_FILE CONFIG_FILE NPOLY MASK_FILE_POLY"
    print "..."
    print "(N) START_W END_W MASK_FILE CONFIG_FILE NPOLY MASK_FILE_POLY"
    print "MIN_DELTA_CHISQ MAX_NITER CUT_MEDIAN_FLUX"
    print "start_w_peak end_w_peak"
    print "wavelength_to_norm width_AA new_back_templates.fits"
    sys.exit(0)
unc_file=sys.argv[1]
clean_file="clean_"+sys.argv[1]
junk_back_list=sys.argv[2]
data=junk_back_list.split(',')
if len(data) == 2:
    back_list=data[0]
    back_list2=data[1]
else:
    back_list=junk_back_list
    back_list2=junk_back_list
outfile=sys.argv[3]
out_elines="elines_"+outfile
out_single="single_"+outfile
out_fit="fit_"+outfile
out_coeffs_file="coeffs_"+outfile
out_fit="output."+outfile+".fits"
out_ps=outfile
#######################################
# Clean previous results
#######################################
call="rm -rf "+outfile
sycall(call)
call="rm -rf "+out_elines
sycall(call)
call="rm -rf "+out_single
sycall(call)
call="rm -rf "+out_fit
sycall(call)
D_SYS_VEL=100
mask_list=sys.argv[4]
config_file=sys.argv[5]
plot=int(sys.argv[6])
if plot == 2:
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.colors as colors
    import matplotlib.cm as cmx
    dev_plot=outfile+".pdf"
    dev_plot_single="single_"+outfile+".pdf"
else:
    if plot == 0:
        matplotlib.use('Agg')
    dev_plot="null"
    dev_plot_single="null"
    import matplotlib.pyplot as plt
    import matplotlib.colors as colors
    import matplotlib.cm as cmx
smooth=1
MIN_CHISQ=1e12
out_file="junk.junk"
factor=1
box=1
deft=0
if len(sys.argv) == 9:
    min=float_(sys.argv[7])
    max=float_(sys.argv[8])
    deft=1
if len(sys.argv) == 11: 
    min=float_(sys.argv[7])
    max=float_(sys.argv[8])
    min_wave=sys.argv[9]
    max_wave=sys.argv[10]
    deft=2
if len(sys.argv) == 12:
    min=float_(sys.argv[7])
    max=float_(sys.argv[8])
    min_wave=sys.argv[9]
    max_wave=sys.argv[10]
    elines_mask=sys.argv[11]
    deft=2
input_redshift=0
if len(sys.argv) == 13:
    min=float_(sys.argv[7])
    max=float_(sys.argv[8])
    min_wave=sys.argv[9]
    max_wave=sys.argv[10]
    elines_mask=sys.argv[11]
    input_redshift=float_(sys.argv[12])
    redshift=input_redshift
    deft=2
f=open(config_file,'r')
line=f.readline()
data=line.split(" ")
data=filter(None,data)
redshift=float_(data[0])
d_redshift=float_(data[1])
min_redshift=float_(data[2])
max_redshift=float_(data[3])
DV=float_(data[4])
RV=float_(data[5])
DS=float_(data[6])
RS=float_(data[7])
MIN_W=float_(data[8])
MAX_W=float_(data[9])
if len(sys.argv) == 16:
    min=float_(sys.argv[7])
    max=float_(sys.argv[8])
    min_wave=sys.argv[9]
    max_wave=sys.argv[10]
    elines_mask=sys.argv[11]
    input_redshift=float_(sys.argv[12])
    input_d_redshift=float_(sys.argv[13])
    input_min_redshift=float_(sys.argv[14])
    input_max_redshift=float_(sys.argv[15])
    redshift=input_redshift
    d_redshift=input_d_redshift
    min_redshift=input_min_redshift
    max_redshift=input_max_redshift
    deft=2
line=f.readline()
data=line.split(" ")
data=filter(None,data)
sigma=data[0]
d_sigma=data[1]
min_sigma=data[2]
max_sigma=data[3]
if len(sys.argv) == 20:
    min=float_(sys.argv[7])
    max=float_(sys.argv[8])
    min_wave=sys.argv[9]
    max_wave=sys.argv[10]
    elines_mask=sys.argv[11]
    input_redshift=float_(sys.argv[12])
    input_d_redshift=float_(sys.argv[13])
    input_min_redshift=float_(sys.argv[14])
    input_max_redshift=float_(sys.argv[15])
    sigma=float_(sys.argv[16])
    d_sigma=float_(sys.argv[17])
    min_sigma=float_(sys.argv[18])
    max_sigma=float_(sys.argv[19])
    redshift=input_redshift
    d_redshift=input_d_redshift
    min_redshift=input_min_redshift
    max_redshift=input_max_redshift
    deft=2
line=f.readline()
data=line.split(" ")
data=filter(None,data)
Av_IN=data[0]
d_Av_IN=data[1]
min_Av=data[2]
max_Av=data[3]
if len(sys.argv) == 24:
    min=float_(sys.argv[7])
    max=float_(sys.argv[8])
    min_wave=sys.argv[9]
    max_wave=sys.argv[10]
    elines_mask=sys.argv[11]
    input_redshift=float_(sys.argv[12])
    input_d_redshift=float_(sys.argv[13])
    input_min_redshift=float_(sys.argv[14])
    input_max_redshift=float_(sys.argv[15])
    sigma=float_(sys.argv[16])
    d_sigma=float_(sys.argv[17])
    min_sigma=float_(sys.argv[18])
    max_sigma=float_(sys.argv[19])
    Av_IN=float_(sys.argv[20])
    d_Av_IN=float_(sys.argv[21])
    min_Av=float_(sys.argv[22])
    max_Av=float_(sys.argv[23])
    redshift=input_redshift
    d_redshift=input_d_redshift
    min_redshift=input_min_redshift
    max_redshift=input_max_redshift
    deft=2
data=min_wave.split(',')
if len(data) == 2:
    min_wave=float_(data[0])
    min_wave2=float_(data[1])
else:
    min_wave=float_(min_wave)
    min_wave2=min_wave
data=max_wave.split(',')
if len(data) == 2:
    max_wave=float_(data[0])
    max_wave2=float_(data[1])
else:
    max_wave=float_(max_wave)
    max_wave2=max_wave
REDSHIFT=redshift
Av_ini=Av_IN
if d_redshift !=0:
    fit_redshift=1
else:
    fit_redshift=0
print "FIT_RED "+str(fit_redshift)+" "+str(d_redshift)+" "+str(len(sys.argv))
line=f.readline()
data=line.split(" ")
data=filter(None,data)
ns=int(data[0])
start_w_min=1e12
end_w_max=-1e12
start_w_E=[]
end_w_E=[]
mask_E=[]
config_E=[]
n_line_E=[]
npoly_E=[]
mask_poly_E=[]
nmin_E=[]
nmax_E=[]
config_line_E=[]
for i in range(0, ns):
    line=f.readline()
    data=line.split(" ")
    data=filter(None,data)
    start_w_e=float_(data[0])
    end_w_e=float_(data[1])
    mask_e=data[2]
    config_e=data[3]
    npoly_e=int(data[4])
    mask_poly_e=data[5]
    nmin_e=float_(data[6])
    nmax_e=float_(data[7])
    start_w_E.extend([start_w_e])
    end_w_E.extend([end_w_e])
    mask_E.extend([mask_e])
    config_E.extend([config_e])
#
# We read all the information
#
    n_line=0
    linef=""
    f2=open(config_e,'r')
    for line in f2:    
        linef=linef+line+";"
        n_line=n_line+1
    config_line_E.extend([linef])
    f2.close()    
    n_line_E.extend([n_line])
    npoly_E.extend([npoly_e])
    mask_poly_E.extend([mask_poly_e])
    nmin_E.extend([nmin_e])
    nmax_E.extend([nmax_e])
    if start_w_e < start_w_min:
        start_w_min=start_w_e
    if end_w_e > end_w_max:
        end_w_max=end_w_e
line=f.readline()
data=line.split(" ")
data=filter(None,data)
MIN_DELTA_CHISQ=float_(data[0])
MAX_NITER=int(data[1])
CUT_MEDIAN_FLUX=float_(data[2])
ABS_MIN=0.5*CUT_MEDIAN_FLUX
line=f.readline()
data=line.split(" ")
data=filter(None,data)
start_w_peak=float_(data[0])
end_w_peak=float_(data[1])
line=f.readline()
data=line.split(" ")
data=filter(None,data)
if len(data) == 3:
    wave_norm=data[0]
    w_wave_norm=data[1]
    new_back_file=data[2]
else:
    wave_norm=[]
    w_wave_norm=[]
    new_back_file=[]
f.close()
#
# SFH template
#
[pdl_flux_c_ini,hdr]=gdata(back_list, 0, header=True)
[nf,n_c]=pdl_flux_c_ini.shape
coeffs=np.zeros([nf,3])
crpix=hdr["CRPIX1"]
cdelt=hdr["CDELT1"]
crval=hdr["CRVAL1"]
n_mc=30
#
# Kinematics template
#
[pdl_flux_c_ini2,hdr2]=gdata(back_list2, 0, header=True)
[nf2,n_c2]=pdl_flux_c_ini2.shape
coeffs2=np.zeros([nf2,3])
crpix2=hdr2["CRPIX1"]
cdelt2=hdr2["CDELT1"]
crval2=hdr2["CRVAL1"]
Av=np.zeros(nf)
d_Av=np.zeros(nf)
for i in range(0, nf):
    Av[i]=Av_IN
    d_Av[i]=d_Av_IN
if mask_list == "none":
    nmask=0
else:
    f=open(mask_list,'r')
    start_mask=[]
    end_mask=[]
    for line in f:
        data=line.split(" ")
        data=filter(None,data)
        if len(data) != 0 and data[0] != "\n":
            start_mask.extend([float_(data[0])])
            end_mask.extend([float_(data[1])])
    nmask=len(start_mask)
    f.close()
    n_mask_org=nmask
if elines_mask == "none":
    nmask_e=0
    nline=0
else:
    f=open(elines_mask,'r')
    nmask_e=0
    nline=0
    w_eline=[]
    start_mask_e=[]
    end_mask_e=[]
    for line in f:
        data=line.split(" ")
        data=filter(None,data)    
        if data[0] != "#":
            w_eline.extend([float_(data[0])])
            start_mask_e.extend([w_eline[nline]*(1+input_redshift)-4*sigma])
            end_mask_e.extend([w_eline[nline]*(1+input_redshift)+4*sigma])
            nmask_e=nmask_e+1
            nline=nline+1
    f.close()
#
#  We read the input spectrum
#
n_unc=0
y_min=1e12
y_max=-1e12
f=open(unc_file,'r')
i_scale=0
FLUX=0
have_error=0
index_unc=[]
wave_unc=[]
flux_unc=[]
flux_unc_org=[]
flux_unc_input=[]
e_flux_unc=[]
color_unc=[]
masked=[]
masked2=[]
masked_Av=[]
flux_masked=[]
flux_masked2=[]
e_flux_unc_kin=[]
wave_scale=0
for line in f:
    data=line.split(' ')
    data=filter(None,data)
    if data[0] != "#":
        index_unc.extend([float_(data[0])])
        wave_unc.extend([float_(data[1])])
        flux_unc.extend([float_(data[2])])
        flux_unc_org.extend([float_(data[2])])
        flux_unc_input.extend([float_(data[2])])
        if len(data) > 2:
        # Variance Column!
            e_flux_unc.extend([np.sqrt(abs(float_(data[3])))]) 
            color_unc.extend([1])#$data[4];
            have_error=1
        else:
            e_flux_unc.extend([np.sqrt(abs(float_(data[2])))/10])
            color_unc.extend([1])
        if np.isnan(flux_unc[n_unc]):
            flux_unc[n_unc]=flux_unc[n_unc-1]
        if flux_unc[n_unc] < y_min:
            y_min=flux_unc[n_unc]
        if flux_unc[n_unc] > y_max:
            y_max=flux_unc[n_unc]
        if n_unc > 0:
            if wave_unc[n_unc-1] <= wave_scale and wave_unc[n_unc] > wave_scale:
                i_scale=n_unc        
        masked.extend([1])
        masked2.extend([1])
        masked_Av.extend([1])
        if flux_unc[n_unc] == 0:
            masked[n_unc]=0    
            masked2[n_unc]=0
        w_test=wave_unc[n_unc-1]
        for j in range(0, nmask):
            if w_test > start_mask[j] and w_test < end_mask[j]:
                masked[n_unc]=0                
                masked2[n_unc]=0        
                masked_Av[n_unc]=0        
        if deft == 2:
            if w_test < min_wave:
                masked[n_unc]=0     
                masked_Av[n_unc]=0
            if w_test > max_wave:
                masked[n_unc]=0       
                masked_Av[n_unc]=0
            if w_test < min_wave2:
                masked2[n_unc]=0       
            if w_test > max_wave2:
                masked2[n_unc]=0        
        for j in range(0, nmask_e):
            if w_test > start_mask_e[j] and w_test < end_mask_e[j]:
                masked2[n_unc]=0     
                masked_Av[n_unc]=0   
        flux_masked.extend([flux_unc[n_unc]*masked[n_unc]])
        flux_masked2.extend([flux_unc[n_unc]*masked2[n_unc]])
        if wave_unc[n_unc] > min_wave and wave_unc[n_unc] < max_wave:
            FLUX=FLUX+flux_masked[n_unc]
        e_flux_unc_kin.extend([e_flux_unc[n_unc]])
        n_unc=n_unc+1
f.close()
sigma_e=np.median(e_flux_unc)
print "SIGMA_E = "+str(sigma_e)
for i in range(0, n_unc):
    if e_flux_unc[i] > 1.5*sigma_e:
        e_flux_unc[i]=1.5*sigma_e    
        e_flux_unc_kin[i]=1.5*sigma_e
if deft == 2:
    y_min=min
    y_max=max
else:
    min_wave=np.amin(wave_unc)
    max_wave=np.amax(wave_unc)
if deft == 1:
    y_min=min
    y_max=max
median_flux=np.median(flux_masked)
dpix_unc=wave_unc[1]-wave_unc[0]
max=3*median_flux
pdl_output=np.zeros([6,n_unc])
#
# We create a kernel
#
med_flux=np.median(flux_unc)
chi_sq_min_now=1e12
min_chi_sq=chi_sq_min_now
ssp_dat=ssp.fit_ssp_lin_no_zero(redshift,sigma,Av,crval2,cdelt2,crpix2,nf2,n_c2,pdl_flux_c_ini2,hdr2,wave_unc,masked_Av,e_flux_unc,flux_unc,n_mc,chi_sq_min_now,min_chi_sq,plot)
min_chi_sq=ssp_dat[0]
print "CUT = "+str(med_flux)+" "+str(CUT_MEDIAN_FLUX)
print str(redshift)+","+str(sigma)
if med_flux < ABS_MIN:
# WHAT TO DO???
# We print all!!!
    sys.exit(0)
if med_flux > CUT_MEDIAN_FLUX:        
    if MIN_W == 0:
        MIN_W = min_wave   
    if MAX_W == 0:
        MAX_W=max_wave
################
# REDSHIFT DETERMINATION
    my_plot=2
    K=0
    nr=0
    chi_r=[]
    red_r=[]
    if d_redshift > 0:
        min_chi_sq=1e30       
        RED=min_redshift
        while RED < max_redshift:
            ssp_dat1=ssp.fit_ssp_lin_no_zero_no_cont(RED,sigma,Av,crval2,cdelt2,crpix2,nf2,n_c2,pdl_flux_c_ini2,hdr2,wave_unc,masked2,e_flux_unc_kin,flux_masked2,n_mc,chi_sq_min_now,min_chi_sq,my_plot)
            chi_now=ssp_dat1[0]
            chi_r.extend([chi_now])
            red_r.extend([RED])
#            print RED,chi_now,d_redshift
            if nr > 1  and chi_r[nr-1] < min_chi_sq and chi_r[nr-1] > 0:
                redshift=red_r[nr-1]    
                min_chi_sq=chi_r[nr-1]    
                K=nr-1
            nr=nr+1
            RED=RED+d_redshift
        #
        # TWO
        #
        e_redshift=d_redshift
        nr=0
        chi_r=[]
        red_r=[]
        RED=redshift-1.5*d_redshift
        max_redshift=redshift+1.5*d_redshift
        d_redshift=0.1*d_redshift
        while RED < max_redshift:
            ssp_dat2=ssp.fit_ssp_lin_no_zero_no_cont(RED,sigma,Av,crval2,cdelt2,crpix2,nf2,n_c2,pdl_flux_c_ini2,hdr2,wave_unc,masked2,e_flux_unc,flux_masked2,n_mc,chi_sq_min_now,min_chi_sq,my_plot)
            chi_now=ssp_dat2[0]
            chi_r.extend([chi_now])
            red_r.extend([RED])
            if nr > 1 and chi_r[nr-1] < chi_r[nr-2] and chi_r[nr-1] < chi_r[nr] and chi_r[nr-1] <= min_chi_sq:
                a=red_r[nr-2]
                b=red_r[nr-1]
                c=red_r[nr]
                fa=chi_r[nr-2]
                fb=chi_r[nr-1]
                fc=chi_r[nr]
                den=(fc-2*fb+fa)
                redshift=c-(b-a)*((fc-fb)/den+0.5)
                slope=abs(0.5*(fc-fb)/(c-b))+abs(0.5*(fa-fb)/(a-b))
                if slope > 0:
                    e_redshift=0.01*redshift/slope
                else:
                    e_redshift=0.01*redshift   
                min_chi_sq=chi_r[nr-1]    
                K=nr-1   
            nr=nr+1
            a_rnd=ran.rand(2)
            RED=RED+d_redshift*(a_rnd[0])
        fit_redshift=0
        d_redshift=0
    else:
        fit_redshift=0
        e_redshift=0
    print "REDSHIFT = "+str(redshift)+" +- "+str(e_redshift)
    #sys.exit()
    REDSHIFT=redshift
    # sigma DETERMINATION    
    K=0
    nr=0
    chi_s=[]
    sigma_s=[]
    print "D_SIGMA = "+str(d_sigma)
    if d_sigma > 0:
        min_chi_sq = 1e30
        SIGMA=min_sigma
        while SIGMA < max_sigma: 
            ssp_dat3=ssp.fit_ssp_lin_no_zero_no_cont(redshift,SIGMA,Av,crval2,cdelt2,crpix2,nf2,n_c2,pdl_flux_c_ini2,hdr2,wave_unc,masked2,e_flux_unc,flux_masked2,n_mc,chi_sq_min_now,min_chi_sq,my_plot)
            chi_now=ssp_dat3[0]
            chi_s.extend([chi_now])
            sigma_s.extend([SIGMA])
            if chi_s[nr-1] < min_chi_sq:
                sigma=sigma_s[nr-1]    
                min_chi_sq=chi_s[nr-1]    
                K=nr
            nr=nr+1
            SIGMA=SIGMA+d_sigma  
        SIGMA=sigma-1.5*d_sigma
        max_sigma=sigma+1.5*d_sigma
        d_sigma=0.33*d_sigma
        #
        #TWO
        #
        nr=0
        chi_s=[]
        sigma_s=[]
        e_sigma=d_sigma
        while SIGMA < max_sigma:
            ssp_dat4=ssp.fit_ssp_lin_no_zero_no_cont(redshift,SIGMA,Av,crval2,cdelt2,crpix2,nf2,n_c2,pdl_flux_c_ini2,hdr2,wave_unc,masked2,e_flux_unc,flux_masked2,n_mc,chi_sq_min_now,min_chi_sq,my_plot)
            chi_now=ssp_dat4[0]
            chi_s.extend([chi_now])
            sigma_s.extend([SIGMA])
            if nr > 1 and chi_s[nr-1] < chi_s[nr-2] and chi_s[nr-1] < chi_s[nr] and chi_s[nr-1] <= min_chi_sq:
                a=sigma_s[nr-2]
                b=sigma_s[nr-1]
                c=sigma_s[nr]
                fa=chi_s[nr-2]
                fb=chi_s[nr-1]
                fc=chi_s[nr]
                den=(fc-2*fb+fa)
                sigma=c-(b-a)*((fc-fb)/den+0.5)
                min_chi_sq=chi_s[nr-1]    
                K=nr
                SIGMA=max_sigma
            nr=nr+1
            a_rnd=ran.rand(2)
            SIGMA=SIGMA+d_sigma*(a_rnd[0])
        slope=(chi_s[nr-1]-min_chi_sq)/(sigma_s[nr-1]-sigma)
        if slope > 0:
            e_sigma=sigma/slope/10.
        else:
            e_sigma=sigma/10.
        fit_sigma=0
        d_sigma=0
    else:
        fit_sigma=0
        e_sigma=0
    sigma=abs(sigma)
    e_sigma=abs(e_sigma)
    print "SIGMA = "+str(sigma)+"+-"+str(e_sigma)
else:
#
# Below the cut!
#
    for i in range(0, nf):
        Av[i]=0
        d_Av[i]=0   
# Av DETERMINATION
K=0
nr=0
chi_Av=[]
Av_s=[]
Av_p_chi=[]
print "D_Av = "+str(d_Av_IN)
nr_min=0
if d_Av_IN > 0:
    min_chi_sq = 1e30
    Av_NOW=min_Av
    while Av_NOW < max_Av:
        for i in range(0, nf):
            Av[i]=Av_NOW
        #
        # Not allow negative coeffs!!!!
        #
        ssp_dat5=ssp.fit_ssp_lin_no_zero(redshift,sigma,Av,crval2,cdelt2,crpix2,nf2,n_c2,pdl_flux_c_ini2,hdr2,wave_unc,masked_Av,e_flux_unc,flux_masked,n_mc,chi_sq_min_now,min_chi_sq,my_plot)
        chi_now=ssp_dat5[0]
        chi_Av.extend([chi_now])
        Av_s.extend([Av_NOW])
        if chi_now > 0:
            Av_p_chi.extend([Av_NOW/(chi_now)])           
        if K == 0 and chi_Av[nr] < min_chi_sq:
            Av_F=Av_s[nr]
            nr_min=nr
            min_chi_sq=chi_now
        if nr > 1 and chi_Av[nr-1] < chi_Av[nr-2] and chi_Av[nr-1] < chi_Av[nr] and chi_Av[nr-1] <= min_chi_sq:
            a=Av_s[nr-2]
            b=Av_s[nr-1]
            c=Av_s[nr]
            fa=chi_Av[nr-2]
            fb=chi_Av[nr-1]
            fc=chi_Av[nr]
            den=(fc-2*fb+fa)
            Av_F=c-(b-a)*((fc-fb)/den+0.5)
            min_chi_sq=chi_Av[nr-1]    
            K=nr
        nr=nr+1
        a_rnd=ran.rand(2)
        Av_NOW=Av_NOW+d_Av_IN*(a_rnd[0])
    if Av_s[nr-1] != Av_F:
        slope=(chi_Av[nr-1]-min_chi_sq)/(Av_s[nr-1]-Av_F)
        if slope > 0 :
            e_Av=abs(Av_F/slope/3.)
        else:
            e_Av=d_Av_IN    
    else:
        e_Av=d_Av_IN
    fit_Av=0
    d_Av_NOW=0
else:
    fit_Av=0     
if d_Av_IN == 0:
    Av_F=Av_IN
if e_Av == 0:
    e_Av=d_Av_IN
print "AV = "+str(Av_F)+" +- "+str(e_Av)
for i in range(0, nf):
    Av[i]=Av_F        
fit_redshift=0
redshift_abs=redshift
delta_chi=10
NITER=0
niter_tmp_max=10
chi_sq_min_now=1e12
min_chi_sq_limit=min_chi_sq
n_mc=10
pdl_rat_master=np.ones(n_unc+1)
[min_chi_sq,pdl_age_mod,pdl_met_mod,pdl_ml,pdl_Av,coeffs,coeffs_N,coeffs_NM,pdl_model_spec_min,pdl_res]=ssp.fit_ssp_lin_no_zero(redshift,sigma,Av,crval2,cdelt2,crpix2,nf2,n_c2,pdl_flux_c_ini2,hdr2,wave_unc,masked_Av,e_flux_unc,flux_unc,n_mc,chi_sq_min_now,min_chi_sq,plot)
#
# We substract the continuum!
# 
pdl_mod_SSP=pdl_model_spec_min
pdl_res_SSP=pdl_res
nx=n_unc
i0_now=int(0.4*nx)
i1_now=int(0.6*nx)
stats_res=np.std(pdl_res[i0_now:i1_now])+np.mean(pdl_res[i0_now:i1_now])
stats_mod=np.mean(pdl_model_spec_min[i0_now:i1_now])
SN=0
if stats_res > 0:
    SN=stats_mod/stats_res
print "Signal-to-Noise = "+str(SN)
old=1       
if old == 1 and SN > 10:
    pdl_model_spec_min[np.where(pdl_model_spec_min == 0)[0]]=1.
    pdl_rat=pdl_res/pdl_model_spec_min+1
    rat=pdl_rat
    med_rat=my.median_filter(int(5*2.354*sigma),rat)
    pdl_med_rat=np.array(med_rat)
    n_unc_1=n_unc-1
    pdl_wave_unc=wave_unc[0]+(wave_unc[1]-wave_unc[0])*np.arange(0,n_unc_1)
    med_rat=my.median_filter(int(7*2.354*sigma),rat)
    med_sigma=int(1.5*sigma)
    if med_sigma < 3:
        med_sigma=3
    med_rat_box=my.median_box(med_sigma,med_rat)
    med_wave_box=my.median_box(med_sigma,wave_unc)
    y_rat = interp1d(med_wave_box, med_rat_box,bounds_error=False,fill_value=0.)(wave_unc)        
    if plot > 0:
        out_ps_now="junk2"
        title="ratio"
        ssp.plot_results_min_max(2,wave_unc,[flux_unc,pdl_model_spec_min,pdl_res,pdl_rat,y_rat],out_ps_now,title,-0.2,1.5)
    i0_now=int(0.4*n_unc)
    i1_now=int(0.6*n_unc)
    stats_rat0=np.mean(y_rat[i0_now:i1_now])
    stats_rat1=np.std(y_rat[i0_now:i1_now])+stats_rat0
    if stats_rat0 > 0 and stats_rat1 > 0.02:
        for i in range(0, n_unc):
            val=y_rat[i]
            if val > 0:
                flux_unc[i]=flux_unc[i]/val
print "Deriving SFH...."
#my pdl_mod_JOINT;
#my pdl_res_JOINT;#=$pdl_res_SSP; #Modificacion 12 de Marzo de 2015 (en caso de no entrar en el ciclo while)
#my pdl_no_gas;
#my age_min;
#my met_min;
#my Av_min;
#my age_min_mass;
#my met_min_mass;
#my Av_min_mass;
coeffs_cat=np.zeros([nf+1,n_mc])
while MIN_CHISQ > MIN_DELTA_CHISQ and NITER < MAX_NITER:
    if NITER == 1:
        MIN_CHISQ=1e12
######################################################################
# Fitting the emission lines
######################################################################
    a_fixed=np.zeros([1,9])
    a_type_fixed=[]
    n_mod_fixed=0
    if ns > 0:
        ks=0
        SYS_VEL=vel_light*REDSHIFT
        REN=[]
        e_REN=[]
        sycall(call)
        for ist in range(0,ns):
            if red_elines > 0:
                SYS_VEL=vel_light*red_elines
                if ist == 0:
                    SYS_VEL_MAX=vel_light*red_elines+D_SYS_VEL
                    SYS_VEL_MIN=vel_light*red_elines-D_SYS_VEL
                else:
                    SYS_VEL_MAX=vel_light*red_elines+D_SYS_VEL
                    SYS_VEL_MIN=vel_light*red_elines-D_SYS_VEL        
            else:
                SYS_VEL=vel_light*REDSHIFT
                if ist == 0:
                    SYS_VEL_MAX=vel_light*REDSHIFT+D_SYS_VEL
                    SYS_VEL_MIN=vel_light*REDSHIFT-D_SYS_VEL
                else:
                    SYS_VEL_MAX=vel_light*REDSHIFT+D_SYS_VEL
                    SYS_VEL_MIN=vel_light*REDSHIFT-D_SYS_VEL
            start_w_e=start_w_E[ist]
            end_w_e=end_w_E[ist]
            mask_e=mask_E[ist]
            config_e=config_E[ist]
            npoly_e=npoly_E[ist]
            mask_poly_e=mask_poly_E[ist]
            nmin_e=nmin_E[ist]
            nmax_e=nmax_E[ist]   
            print "CONF="+config_e
            wave_elines=[]
            flux_elines=[]
            flux_e_elines=[]
            masked_elines=[]
            n_e=0
            for i in range(0, n_unc):
                if wave_unc[i] > start_w_e and wave_unc[i] < end_w_e:
                    wave_elines.extend([wave_unc[i]])
                    flux_elines.extend([flux_unc_org[i]-pdl_mod_SSP[i]])
                    flux_e_elines.extend([abs(e_flux_unc[i])])
                    masked_elines.extend([1])
                    n_e=n_e+1
            pdl_wave_elines=np.array(wave_elines)
            pdl_flux_elines=np.array(flux_elines)
            pdl_flux_e_elines=np.array(flux_e_elines)
            pdl_masked_elines=np.array(masked_elines)
            stats0=np.mean(pdl_flux_elines)
            stats4=np.amax(pdl_flux_elines)
            y_max=stats4-stats0
            deft=1
            data=filter(None, config_line_E[ist].split(';')[0].split(" "))
            #print float_(filter(None, config_line_E[0].split(';')[4].split(" ")))[2]
            junk=data[0]
            n_mod=int(data[1])
            chi_goal=float_(data[2])
            d_chi_goal=float_(data[3])
            n_line=n_line_E[ist]
            i_mod=1
            typef=[]
            a=np.zeros([n_mod,9])
            ia=np.zeros([n_mod,9])
            ea=np.zeros([n_mod,9])
            a0=np.zeros([n_mod,9])
            a1=np.zeros([n_mod,9])
            link=np.zeros([n_mod,9])
            for ii in range(0, n_mod):
                cnf=filter(None, config_line_E[ist].split(';')[i_mod].split(" "))  
                i_mod=i_mod+1
                typef.extend(cnf)
                for j in range(0, 9):
                    data=config_line_E[ist].split(';')[i_mod].replace('\t',' ')
                    data=filter(None, data.split(' '))
                    i_mod=i_mod+1
                    a[ii][j]=float_(data[0])
                    ia[ii][j]=float_(data[1])
                    ea[ii][j]=0
                    a0[ii][j]=float_(data[2])
                    a1[ii][j]=float_(data[3])
                    link[ii][j]=float_(data[4])
                if deft == 1:
                    a1_max=2*y_max*(a[ii][2]*((2*3.1416)**0.5))
                    a0_min=0.01*1.2*y_max*(a[ii][2]*((2*3.1416)**0.5))
                    if a1[ii][1] > a1_max:
                        a1[ii][1]=a1_max
                        a0[ii][1]=a0_min
            #
            # Force vicitiny in the velocity
            #
            a[0][3]=SYS_VEL
            ia[0][3]=1
            a0[0][3]=SYS_VEL_MIN    
            a1[0][3]=SYS_VEL_MAX
            i_ter=0
            n_loops=5
            n_mc_e=30
            breakt=0
            scale_ini=0.15;
            deft=0;
            pdl_model=np.zeros(n_e)
            pdl_model_cont=np.zeros(n_e)
            pdl_model_tmp=np.zeros(n_e)
            pdl_model_cont_tmp=np.zeros(n_e)
            a_out_now=ssp.copy_a(n_mod,a)
            a_out_tmp=ssp.copy_a(n_mod,a)
            chi_sq_now=1e12
            a_results=np.zeros([1, n_mod, 9])
            nnk=0
            max_time=5
            time=ssp.get_seconds()
            d_time=ssp.get_seconds()-time
            i_loops=0
            ran.seed(None)
            while i_ter < n_loops and breakt == 0:
                chi_iter=chi_sq_now
                chi_single=0
                [chi_sq_now,pdl_a,pdl_model_tmp,pdl_model_cont_tmp]=ssp.fit_elines_grad_rnd_new(pdl_wave_elines,pdl_flux_elines,pdl_flux_e_elines,n_mod,chi_goal,d_chi_goal,typef,a_out_tmp,ia,a0,a1,link,n_mc_e,pdl_masked_elines,deft,scale_ini)#,max_time)
                a_out_now=ssp.copy_a_pdl(n_mod,pdl_a)
                #print chi_sq_now, pdl_a[:,1],a_out_tmp[:,1]
                if chi_sq_now < chi_iter:
#####################################################                    
                # Close to a result, narrow the range                                                                              
                    for i in range(0, n_mod):
                        for j in range(0, 9):
                            if typef[i] == "eline\n":
                                if ia[i][j] == 1:
                                    if link[i][j] == -1:
                                        delta_now=abs(a1[i][j]-a0[i][j])/(2.)
                                        a0_tmp=a0[i][j]
                                        a1_tmp=a1[i][j]
                                        if j != 3:
                                            a0_tmp=a_out_now[i][j]-delta_now
                                            a1_tmp=a_out_now[i][j]+delta_now
                                        else:
                                            a0_tmp=a_out_now[i][j]-0.5*delta_now
                                            a1_tmp=a_out_now[i][j]+0.5*delta_now
                                        if a0_tmp < a0[i][j]:
                                            a0_tmp=a0[i][j]
                                        if a1_tmp > a1[i][j]:
                                            a1_tmp=a1[i][j]
                                        a0[i][j]=a0_tmp
                                        a1[i][j]=a1_tmp
####################################################                             
                    a_out_tmp=ssp.copy_a(n_mod,a_out_now)
                    a_results=ssp.copy_a_results(n_mod,nnk,a_out_now,a_results)    
                    pdl_model=pdl_model_tmp
                    pdl_model_cont=pdl_model_cont_tmp
                    nnk=nnk+1
                    i_ter=i_ter+1
                else:
                    rnd_a=ran.rand(10);
                    a_out_now=ssp.copy_a(n_mod,a_out_now)
                    i_loops=i_loops+1
                    if i_loops > 5*n_loops:
                        breakt=1
            out_ps_now="fit_"+outfile+"."+str(start_w_e)+"_"+str(end_w_e)
            title="["+str(start_w_e)+","+str(end_w_e)+"]"
            if pdl_model.shape[0] == len(pdl_wave_elines):
                pdl_model=np.transpose(pdl_model)
            ssp.plot_results(plot,pdl_wave_elines,[pdl_flux_elines,pdl_model[0,:],(pdl_flux_elines-pdl_model[0,:])],out_ps_now,title)
            print "----------------------------------------";
            a_final=ssp.mean_a_results_last(n_mod,nnk,a_results,ia)
            #
            # Background noise                                                                                                         
            #
            pdl_res_now=pdl_flux_elines-pdl_model
            stats_back1=np.mean(pdl_res_now)+np.std(pdl_res_now)
            a_final=ssp.add_back_noise(n_mod,a_final,typef,chi_sq_now,stats_back1)
            ssp.print_a_final(n_mod,a_final,typef,chi_sq_now)
            out_fit_spectra=out_elines
            ssp.print_a_final_file_add(n_mod,a_final,typef,chi_sq_now,out_fit_spectra)  
            [n_mod_fixed,junk_a_fixed,junk_a_type_fixed]=ssp.add_a_results_elines(n_mod,a_final,typef,n_mod_fixed,a_fixed,a_type_fixed)
            a_fixed=junk_a_fixed
            a_type_fixed=junk_a_type_fixed
            nmin_e=int(0.1*n_unc)
            nmax_e=int(0.9*n_unc)
###############################
            # Low order polynomical!
            out_fit_now=out_fit+"."+str(start_w_e)+"_"+str(end_w_e)+".pdf"
            box=int(sigma*6)
            print "DONE FIT ELINES CONFIG "+str(ist)
        #
        # We create a FIXED model of the emission lines
        #
        pdl_model_elines=np.zeros(n_unc)
        pdl_model_cont=np.zeros(n_unc)
        pdl_wave_elines=np.array(wave_unc)
        NN=len(pdl_wave_elines)
        NN1=len(pdl_model_elines)
        for i in range(0, n_mod_fixed):
            pdl_tmp=ssp.create_single_model(pdl_wave_elines,i,a_type_fixed,a_fixed)
            NN2=len(pdl_tmp[0,:])
            pdl_model_elines=pdl_model_elines+pdl_tmp[0,:]
        #
        # We remove the gas before a new iteration
        #
        for i in range(0, n_unc):
            flux_unc[i]=flux_unc_org[i]-pdl_model_elines[i]
        pdl_mod_JOINT=pdl_mod_SSP+pdl_model_elines
        pdl_res_JOINT=pdl_res_SSP-pdl_model_elines
        pdl_no_gas=np.array(flux_unc)
        #############################################################
        # We rescale!
        ##############################################################
        y_rat=np.ones(nx+1)
        jy=0
        if SN > 10:
            pdl_mod_JOINT[np.where(pdl_mod_JOINT == 0)[0]]=1.
            pdl_rat=pdl_res_JOINT/pdl_mod_JOINT+1
            rat=pdl_rat
            n_unc_1=n_unc-1
            pdl_wave_unc=wave_unc[0]+(wave_unc[1]-wave_unc[0])*np.arange(0, n_unc_1)
            med_rat=my.median_filter(int(5*2.354*sigma),rat);
            med_sigma=int(1.5*sigma)
            if med_sigma < 3:
                med_sigma=3
            med_rat_box=my.median_box(med_sigma,med_rat)
            med_wave_box=my.median_box(med_sigma,wave_unc) 
            y_rat = interp1d(med_wave_box, med_rat_box,bounds_error=False,fill_value=0.)(wave_unc) 
            i0_now=int(0.4*nx)
            i1_now=int(0.6*nx)
            stats_rat0=np.mean(y_rat[i0_now:i1_now])
            stats_rat1=np.std(y_rat[i0_now:i1_now])+stats_rat0        
            if plot > 1:
                out_ps_now="junk3"
                title="ratio = "+str(stats_rat0)+", rms="+str(stats_rat1)
                print title
                ssp.plot_results_min_max(2,wave_unc,[flux_unc,pdl_model_spec_min,pdl_res,pdl_rat,y_rat],out_ps_now,title,-0.2,1.5)
            if stats_rat0 > 0 and stats_rat1 > 0.02:
                if jy == 0:
                    # Continuum shape correction on/off
                    pdl_rat_master=y_rat
                    pdl_rat_master[np.where(pdl_rat_master == 0)[0]]=1.
                    y_rat=pdl_rat_master
                else:
                    y_rat=pdl_rat_master
        for i in range(0, n_unc):
            val=y_rat[i]
            if val > 0:
                flux_unc[i]=flux_unc[i]/val
                flux_unc_org[i]=flux_unc_org[i]/val
        ##############################################################
        # End re-scale
        ##############################################################    
        ML=0
        if med_flux > CUT_MEDIAN_FLUX:
            n_mc=20
            [min_chi_sq,pdl_age_mod,pdl_met_mod,pdl_ml,pdl_Av,coeffs,coeffs_N,coeffs_NM,pdl_mod_SSP,pdl_res_SSP,coeffs_N_input,e_coeffs_N_input]=ssp.fit_ssp_lin_MC(redshift,sigma,Av,crval,cdelt,crpix,nf,n_c,pdl_flux_c_ini,hdr,wave_unc,masked,e_flux_unc,flux_unc,n_mc,chi_sq_min_now,MIN_CHISQ,plot)
            smooth_ratiot=ssp.smooth_ratio(flux_unc,pdl_mod_SSP,int(sigma))
            pdl_mod_SSP_no_cor=np.copy(pdl_mod_SSP)
            pdl_mod_SSP=pdl_mod_SSP*smooth_ratiot
            f1=open(out_coeffs_file, "w")
            f1.write("# ID   AGE     MET    COEFF   Min.Coeff  log(M/L)   AV   N.Coeff   Err.Coeff\n")
            print "------------------------------------------------------------------------------"
            print "ID   AGE     MET    COEFF   Min.Coeff  log(M/L)   AV   N.Coeff   Err.Coeff"
            print "------------------------------------------------------------------------------"
            age_mod=pdl_age_mod
            met_mod=pdl_met_mod
            Av_mod=pdl_Av
            ml=pdl_ml
            a_coeffs=coeffs[:,0]
            a_e_coeffs=coeffs[:,1]
            a_min_coeffs=coeffs[:,2]
            a_coeffs_N=coeffs_N  
            a_e_coeffs_N=a_e_coeffs
            l_age_min=0
            l_met_min=0
            l_Av_min=0
            l_age_min_mass=0
            l_met_min_mass=0
            l_Av_min_mass=0
            e_l_age_min=0
            e_l_met_min=0
            e_l_Av_min=0
            e_l_age_min_mass=0
            e_l_met_min_mass=0
            e_l_Av_min_mass=0
            for k in range(0, nf):
                if a_coeffs[k] > 0:
                    a_e_coeffs_N[k]=a_e_coeffs[k]*(a_coeffs_N[k]/a_coeffs[k])
                else:
                    a_e_coeffs_N[k]=0
                l_age_min=l_age_min+a_coeffs[k]*np.log10(age_mod[k])
                l_met_min=l_met_min+a_coeffs[k]*np.log10(met_mod[k])
                l_Av_min=l_Av_min+a_coeffs[k]*np.log10(Av_mod[k])                
                l_age_min_mass=l_age_min_mass+ml[k]*a_coeffs_N[k]*np.log10(age_mod[k])
                l_met_min_mass=l_met_min_mass+ml[k]*a_coeffs_N[k]*np.log10(met_mod[k])
                l_Av_min_mass=l_Av_min_mass+ml[k]*a_coeffs_N[k]*np.log10(Av_mod[k])
                e_l_age_min=e_l_age_min+a_e_coeffs[k]*np.log10(age_mod[k])
                e_l_met_min=e_l_met_min+a_e_coeffs[k]*np.log10(met_mod[k])
                e_l_Av_min=e_l_Av_min+a_e_coeffs[k]*np.log10(Av_mod[k])              
                e_l_age_min_mass=e_l_age_min_mass+ml[k]*a_e_coeffs_N[k]*np.log10(age_mod[k])
                e_l_met_min_mass=e_l_met_min_mass+ml[k]*a_e_coeffs_N[k]*np.log10(met_mod[k])
                e_l_Av_min_mass=e_l_Av_min_mass+ml[k]*a_e_coeffs_N[k]*np.log10(Av_mod[k])
                ML=ML+ml[k]*a_coeffs_N[k]
                C_ini=coeffs_N_input[k]
                e_C_ini=e_coeffs_N_input[k]                  
                f1.write(("%2d" % k)+" "+("%7.4f" % age_mod[k])+" "+("%7.4f" % met_mod[k])+" "+("%7.4f" % a_coeffs_N[k])+" "+("%7.4f" % a_min_coeffs[k])+" "+("%4.4f" % np.log10(ml[k]))+" "+("%4.2f" % Av_mod[k])+" "+("%7.4f" % a_coeffs[k])+" "+("%7.4f" % a_e_coeffs[k])+"\n")
                if a_coeffs[k] > 1e-5:
                    print ("%2d" % k)+" "+("%7.4f" % age_mod[k])+" "+("%7.4f" % met_mod[k])+" "+("%7.4f" % a_coeffs_N[k])+" "+("%7.4f" % a_min_coeffs[k])+" "+("%4.4f" % np.log10(ml[k]))+" "+("%4.2f" % Av_mod[k])+" "+("%7.4f" % a_coeffs[k])+" "+("%7.4f" % a_e_coeffs[k])+" "+("%7.4f" % C_ini)+" "+("%7.4f" % e_C_ini)
            print "------------------------------------------------------------------------------"
            f1.close
            age_min=10**(l_age_min)
            met_min=10**(l_met_min)
            Av_min=10**(l_Av_min)
            age_min_mass=10**(l_age_min_mass/ML)
            met_min_mass=10**(l_met_min_mass/ML)
            Av_min_mass=10**(l_Av_min_mass/ML)
            e_age_min=abs(0.43*e_l_age_min*age_min)
            e_met_min=abs(0.43*e_l_met_min*met_min)
            e_Av_min=abs(0.43*e_l_Av_min*Av_min)
            e_age_min_mass=abs(0.43*e_l_age_min*age_min_mass)
            e_met_min_mass=abs(0.43*e_l_met_min*met_min_mass)
            e_Av_min_mass=abs(0.43*e_l_Av_min*Av_min_mass)
            if min_chi_sq > 0:
                delta_chi=abs((chi_sq_min_now-min_chi_sq)/min_chi_sq)
        wpeak=6562
        Fpeak=-1e12
        pdl_mod_JOINT=pdl_mod_SSP+pdl_model_elines
        pdl_res_JOINT=pdl_res_SSP-pdl_model_elines
        pdl_no_gas=np.array(flux_unc)
        # Copy output!
        pdl_output[0,:]=np.array(flux_unc_org)
        pdl_output[1,:]=pdl_mod_SSP
        pdl_output[2,:]=pdl_mod_JOINT
        pdl_res_SSP=np.array(flux_unc_org)-pdl_mod_SSP
        pdl_res_SSP_no_cor=np.array(flux_unc_input)-pdl_mod_SSP_no_cor
        pdl_output[3,:]=pdl_res_SSP_no_cor
        pdl_tmp=np.array(flux_unc_org)
        nx_1=n_unc#-1
        if len(pdl_rat_master)-len(pdl_mod_JOINT)==1:
            pdl_res_JOINT=pdl_tmp/(pdl_rat_master[0:nx_1])-pdl_mod_JOINT
        else:
            pdl_res_JOINT=pdl_tmp/(pdl_rat_master)-pdl_mod_JOINT
        pdl_output[4,:]=pdl_res_JOINT
        pdl_output[5,:]=np.array(flux_unc_org)-(pdl_res_SSP-pdl_res_JOINT)
        title="X="+str(chi_sq_now)+" T="+str(age_min)+" ("+str(age_min_mass)+") Z="+str(met_min)+" ("+str(met_min_mass)+") Av="+str(Av_min)+" z="+str(redshift)+" sigma="+str(sigma)
        ssp.plot_results(plot,pdl_wave_elines,pdl_output,out_ps,title)
    print "I.Iter = "+str(NITER)+" DONE"
    NITER=NITER+1
# Write output file
h=pyf.PrimaryHDU().header
h["NAXIS"]=2 
h["NAXIS1"]=n_unc 
h["NAXIS2"]=6
h["COMMENT"]="OUTPUT auto_ssp_elines_rnd.pl FITs"
h["CRVAL1"]=wave_unc[0]
h["CDELT1"]=wave_unc[1]-wave_unc[0];
h["CRPIX1"]=1
if pt.exists(out_fit) == False:
    wfits(out_fit,pdl_output,h)
else:
    sycall("rm "+out_fit)
    wfits(out_fit,pdl_output,h)
################################
print "--------------------------------------------------------------"
pdl_masked=np.array(masked)
pdl_chi_now=((pdl_masked*pdl_res_JOINT)**2)/((np.array(e_flux_unc))**2)
pdl_chi_now[np.isnan(pdl_chi_now)]=0
chi_joint=np.sum(pdl_chi_now)
chi_joint=(chi_joint/(n_unc-n_mod_fixed-nf-1))**0.5
rms=np.std(pdl_masked*pdl_res_JOINT)
j1=int(0.4*n_unc)
j2=int(0.6*n_unc)
rms=np.std(pdl_res_JOINT[j1:j2])
pdl_flux_unc_now=np.array(flux_unc)
med_flux=np.median(pdl_flux_unc_now[j1:j2])
title="X="+str(chi_joint)+" T="+str(age_min)+" ("+str(age_min_mass)+") Z="+str(met_min)+" ("+str(met_min_mass)+") Av="+str(Av_min)+" z="+str(redshift)+" sigma="+str(sigma)
ssp.plot_results(plot,wave_unc,pdl_output,out_ps,title)
MASS=ML*med_flux
lML=np.log10(ML)
print "MSP CHISQ="+str(chi_joint)+" AGE="+str(age_min)+"+-"+str(e_age_min)+" MET="+str(met_min)+"+-"+str(e_met_min)+" AV="+str(Av_min)+"+-"+str(e_Av_min)+" REDSHIFT="+str(redshift)+"+-"+str(e_redshift)+" SIGMA_DISP="+str(sigma)+"+-"+str(e_sigma)+" RMS="+str(rms)+" MED_FLUX="+str(med_flux)+" AGE_mass="+str(age_min_mass)+"+-"+str(e_age_min_mass)+" MET_mass="+str(met_min_mass)+"+-"+str(e_met_min_mass)+" MASS="+str(MASS)+" log_M/L="+str(lML)
j1=int(0.4*n_unc)
j2=int(0.6*n_unc)
wave_norm=(wave_unc[j1]+wave_unc[j2])/2.
f=open(outfile, "w")
f.write("# (1) MIN_CHISQ\n")
f.write("# (2) LW Age (Gyr)\n")
f.write("# (3) LW Age error\n")
f.write("# (4) LW metallicity\n")
f.write("# (5) LW metallicity error\n")
f.write("# (6) Av\n")
f.write("# (7) AV error\n")
f.write("# (8) redshift \n")
f.write("# (9) redshift error\n")
f.write("# (10) velocity dispersion sigma, in AA\n")
f.write("# (11) velocity dispersion error\n")
f.write("# (12) median_FLUX\n")
f.write("# (13) redshift_ssp\n")
f.write("# (14) med_flux \n")
f.write("# (15) StdDev_residual \n")
f.write("# (16) MW Age (Gyr)\n")
f.write("# (17) MW Age error\n")
f.write("# (18) MW metallicity\n")
f.write("# (19) MW metallicity error\n")
f.write("# (20) Systemic Velocity km/s \n")
f.write("# (21) Log10 Average Mass-to-Light Ratio \n")
f.write("# SSP_SFH $back_list \n")
f.write("# SSP_KIN $back_list2 \n")
f.write("# WAVE_NORM $wave_norm AA\n")
if chi_joint == 0:
    chi_joint=1
f.write(str(chi_joint)+","+str(age_min)+","+str(e_age_min)+","+str(met_min)+","+str(e_met_min)+","+str(Av_min)+","+str(e_Av_min)+","+str(redshift)+","+str(e_redshift)+","+str(sigma)+","+str(e_sigma)+","+str(FLUX)+","+str(redshift_abs)+","+str(med_flux)+","+str(rms)+","+str(age_min_mass)+","+str(e_age_min_mass)+","+str(met_min_mass)+","+str(e_met_min_mass)+","+str(SYS_VEL)+","+str(lML)+"\n")
f.close
sec_end=ssp.print_time()
sec_total=sec_end-sec_ini
print "# SECONDS = "+str(sec_total)
#
# Write the output!
# 
#
    