#!/usr/bin/python
from mpi4py import MPI
import sys
import numpy as np
from numpy import float_
from numpy import absolute as abs
from numpy import random as ran
import time
import my_auto_ssp_elines_rnd as ssp
import my_auto_ssp_elines_rnd_sigma_inst as ssp1
import pyfits as pyf
from pyfits import getheader as ghead
from pyfits import getdata as gdata
from pyfits import writeto as wfits
from scipy.interpolate.interpolate import interp1d
import my as my
import os.path as pt
import matplotlib
matplotlib.use('Agg')

comm = MPI.COMM_WORLD
nprocs = comm.Get_size()
rank   = comm.Get_rank()

def sycall(comand,si=0):
    from subprocess import call
    line=comand.split(" ")
    fcomand=[]
    fcomand.extend(line)
    linp=''
    nx=len(fcomand)
    for i in range(0, nx):
        linp=linp+fcomand[i]+" "
    if si == 0:
        print linp
    call(fcomand)
    
sys.argv=filter(None,sys.argv)
ran.seed(None)
vel_light=299792.458
red_elines=0.0
if rank == 0:
    sec_ini=ssp.print_time()
    time1=ssp.get_time()
    if len(sys.argv) < 7:
        print "USE: auto_ssp_elines_rnd_rss_sigma_inst_guided_mpi.py OUTFIT3D_SPEC1.RSS.fits,ERROR.RSS.fits SSP_SFH.fits,INPUT_KIN.csv,SIGMA_INST OUTFILE MASK_LIST CONFIG_FILE PLOT [min max] [wmin wmax] [redshift_elines_to_mask] [input_redshift delta_redshift min_redshift max_redshift] [input_sigma delta_sigma min_sigma max_sigma] [input_Av delta_Av min_Av max_Av]"
        print "CONFIG_FILE:"
        print "redshift delta_redshift min_redshift max_redshift"
        print "sigma delta_sigma min_sigma max_sigma (km/h)"
        print "Av delta_Av min_Av max_Av [Same range for all - magnitude ]"
        print "N_SYSTEMS"
        print "(1) START_W END_W MASK_FILE CONFIG_FILE NPOLY MASK_FILE_POLY"
        print "..."
        print "(N) START_W END_W MASK_FILE CONFIG_FILE NPOLY MASK_FILE_POLY"
        print "MIN_DELTA_CHISQ MAX_NITER CUT_MEDIAN_FLUX"
        print "start_w_peak end_w_peak"
        print "wavelength_to_norm width_AA new_back_templates.fits"
        sys.exit(0)
else:
    time1=ssp.get_time()
    if len(sys.argv) < 7:
        sys.exit(0)
n_mc=30
JUNK_input_rss=sys.argv[1]
data=JUNK_input_rss.split(',')
if len(data) == 2:
    input_rss=data[0]
    input_e_rss=data[1]
else:
    input_rss=JUNK_input_rss
    input_e_rss="none"
clean_file="clean_"+sys.argv[1]
junk_back_list=sys.argv[2]
data=junk_back_list.split(',')
if len(data) == 4:
    back_list=data[0]
    back_list_kin=data[1]
    back_list2=data[2]
    sigma_inst=float_(data[3])
else:
    back_list=junk_back_list
    back_list2=junk_back_list
    sigma_inst=0.01
if sigma_inst == 0:
    sigma_inst=0.01
outfile=sys.argv[3]
out_elines="elines_"+outfile
out_single="single_"+outfile
out_fit="fit_"+outfile
out_coeffs_file="coeffs_"+outfile
out_fit="output."+outfile+".fits"
out_ps="pdf_"+outfile
#######################################
# Clean previous results
#######################################
if rank == 0:
    call="rm -rf "+outfile
    sycall(call,1)
    call="rm -rf "+out_elines
    sycall(call,1)
    call="rm -rf "+out_single
    sycall(call,1)
    call="rm -rf "+out_fit
    sycall(call,1)
D_SYS_VEL=300
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
sigma=float_(data[0])
d_sigma=float_(data[1])
min_sigma=float_(data[2])
max_sigma=float_(data[3])
sigma_mean=np.sqrt(sigma_inst**2.0+(5000.0*(sigma/300000.0))**2.0)
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
if rank == 0:
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
crpix1=hdr["CRPIX1"]
cdelt1=hdr["CDELT1"]
crval1=hdr["CRVAL1"]
#
# Kinematics template
#

[pdl_flux_c_ini_kin, hdr2]=gdata(back_list_kin, 0, header=True)
[nf_kin,n_c_kin]=pdl_flux_c_ini_kin.shape
coeffs_kin=np.zeros([nf_kin,3])
crpix_kin=hdr2["CRPIX1"]
cdelt_kin=hdr2["CDELT1"]
crval_kin=hdr2["CRVAL1"]
nf2=0
Av_fix=[]
red_fix=[]
sigma_fix=[]
e_Av_fix=[]
e_red_fix=[]
e_sigma_fix=[]
#print back_list2
#sys.exit()
fin=open(back_list2,"r")
for line in fin:
    if not "#" in line:
        data=line.replace("\n","").split(",")
        data=filter(None,data)
        Av_fix.extend([float_(data[5])])
        red_fix.extend([float_(data[7])])
        sigma_fix.extend([float_(data[9])])
        e_Av_fix.extend([float_(data[6])])
        e_red_fix.extend([float_(data[8])])
        e_sigma_fix.extend([float_(data[10])])
        nf2=nf2+1
fin.close()
Av_fix=np.array(Av_fix)
red_fix=np.array(red_fix)
sigma_fix=np.array(sigma_fix)
e_Av_fix=np.array(e_Av_fix)
e_red_fix=np.array(e_red_fix)
e_sigma_fix=np.array(e_sigma_fix)
print back_list2+" readed, "+str(nf2)+" values found"
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
            start_mask_e.extend([w_eline[nline]*(1+input_redshift)-4*sigma_mean])
            end_mask_e.extend([w_eline[nline]*(1+input_redshift)+4*sigma_mean])
            nmask_e=nmask_e+1
            nline=nline+1
    f.close()
#
#  We read the input spectrum
#
[pdl_input_cube, hd]=gdata(input_rss, 0, header=True)
nx=hd["NAXIS2"]
NY=hd["NAXIS3"]
crpix=hd["CRPIX1"]
crval=hd["CRVAL1"]
cdelt=hd["CDELT1"]
pdl_input_rss=pdl_input_cube[:,:,5]
if input_e_rss != "none":
    pdl_input_e_rss=gdata(input_e_rss, 0)
    pdl_input_e_rss[np.isnan(pdl_input_e_rss)]=1
pdl_wave_org=crval+cdelt*(np.arange(0, nx)+1-crpix)
wave_unc=pdl_wave_org
nm=nx
pdl_masked=np.ones(nm)
org=np.zeros([NY,nm])
mod=np.zeros([NY,nm])
res=np.zeros([NY,nm])
masked=np.zeros(nx)
masked2=np.zeros(nx)
masked_Av=np.zeros(nx)
for i in range(0, nx):
    w_test=pdl_wave_org[i]
    masked[i]=1
    masked2[i]=1
    masked_Av[i]=1
    for j in range(0, nmask):
        if w_test > start_mask[j] and w_test < end_mask[j]:
            masked[i]=0        
            masked2[i]=0        
            masked_Av[i]=0
        if deft == 2:
            if w_test < min_wave:
                masked[i]=0        
                masked_Av[i]=0
            if w_test > max_wave:
                masked[i]=0     
                masked_Av[i]=0   
            if w_test < min_wave2:
                masked2[i]=0
            if w_test > max_wave2:
                masked2[i]=0
    for j in range(0, nmask_e):
        if w_test > start_mask_e[j] and w_test < end_mask_e[j]:
            masked2[i]=0     
            masked_Av[i]=0   
pdl_masked_org=masked
pdl_masked2_org=masked2
n_unc=nx
if deft == 2:
    y_min=min
    y_max=max
else:
    min_wave=np.amin(wave_unc)
    max_wave=np.amax(wave_unc)
if deft == 1:
    y_min=min
    y_max=max
#
# We check the median value  
#
#median_flux=np.median(flux_masked)
dpix_unc=wave_unc[1]-wave_unc[0]
#max=3*median_flux
pdl_output=np.zeros([NY,n_unc,6])
pdl_output_now=np.zeros([n_unc,6])
if rank == 0:
    f=open(outfile, "a")
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
    f.write("# (21) log10 Average Mass-to-Light Ratio\n")
    f.write("# (22) row of the RSS\n")
    f.write("# SSP_SOFH "+back_list+" \n")
    f.write("# SSP_KIN "+back_list2+" \n")
    f.write("# WAVE_NORM "+str(wave_norm)+" AA\n")
    f.write("# INPUT_RSS "+input_rss+"\n")
    f.close
    f=open(out_coeffs_file, "w")
    f.write("# ID   AGE     MET    COEFF   Min.Coeff  log(M/L)   AV   N.Coeff   Err.Coeff\n")
    f.close
INPUT_SIGMA=sigma
INPUT_D_SIGMA=d_sigma
INPUT_MIN_SIGMA=min_sigma
INPUT_MAX_SIGMA=max_sigma
INPUT_REDSHIFT=redshift
INPUT_D_REDSHIFT=d_redshift
INPUT_MIN_REDSHIFT=min_redshift
INPUT_MAX_REDSHIFT=max_redshift
INPUT_AV=Av_IN
age_min=0
met_min=0
Av_min=0
age_min_mass=0
met_min_mass=0
Av_min_mass=0
pdl_rat_master=np.ones(nx+1)
#DISTRIBUTE THE SPECTRAS
e_line_jy=[]
out_coeffs_f=[]
openfile_f=[]
nt=int(NY)/int(nprocs)
nst=np.mod(NY,nprocs)
ds=1
i_in_t=0
i_fi_t=0
if rank == 0:
    ind1=np.zeros(nprocs)
    ind2=np.zeros(nprocs)
for i in range(0, nprocs):
    if nst == i:
        ds=0
    i_in_t=i_fi_t
    i_fi_t=nt+ds+i_in_t
    if rank == 0:
        ind1[i]=i_in_t
        ind2[i]=i_fi_t
    if rank == i:
        i_in=i_in_t
        i_fi=i_fi_t
#BEGIN OF THE RSS LOOP
for jy in range(i_in, i_fi):
    
    sigma=sigma_fix[jy]
    redshift=red_fix[jy]
    Av_IN=Av_fix[jy]
    Av_min=Av_IN
    e_sigma=e_sigma_fix[jy]
    e_redshift=e_red_fix[jy]
    e_Av=e_Av_fix[jy]
    for i in range(0, nf):
        Av[i]=Av_IN
        d_Av[i]=d_Av_IN
    #print Av_fix[jy]
    #sys.exit()
    pdl_flux=pdl_input_rss[jy,:]
    nx_1=nx
    pdl_rat_master[0:nx_1][np.where(pdl_rat_master[0:nx_1] == 0)[0]]=1
    pdl_flux=pdl_flux/(pdl_rat_master[0:nx_1])
    pdl_flux[np.isnan(pdl_flux)]=0
    statst=my.stats(pdl_flux)    
#
#  This has to change!
#
    pdl_e_flux=[]
    if input_e_rss != "none":
        pdl_e_flux=np.sqrt(abs(pdl_input_e_rss[jy,:]))
        statst_p=my.stats(pdl_e_flux)
        pdl_e_flux[np.where(pdl_e_flux <= 0)[0]]=3*statst_p[2]
        pdl_e_flux[np.where(pdl_e_flux >= 3*statst_p[2])[0]]=3*statst_p[2]
    else:
        pdl_e_flux=sqrt(abs(pdl_flux+0.5*statst[1]))/10.
    e_flux_unc=np.array(pdl_e_flux)
    flux_unc=np.array(pdl_flux)
    flux_unc_org=np.array(pdl_flux)
    flux_unc_input=np.array(pdl_flux)
    flux_masked=np.array(pdl_flux*pdl_masked_org)
    flux_masked2=np.array(pdl_flux*pdl_masked2_org)
    med_flux=statst[2]
    chi_sq_min_now=1e12
    #pdl_mod_SSP
    #pdl_res_SSP
    min_chi_sq=1e12
    #ssp_dat0=ssp.fit_ssp_lin_no_zero_no_cont(redshift,sigma,Av,crval2,cdelt2,crpix2,nf2,n_c2,pdl_flux_c_ini2,hdr2,wave_unc,masked_Av,e_flux_unc,flux_unc,n_mc,chi_sq_min_now,min_chi_sq,plot)
    #min_chi_sq=ssp_dat0[0]
    #print "CUT = "+str(med_flux)+" "+str(CUT_MEDIAN_FLUX)
    #print str(redshift)+","+str(sigma)
    if med_flux > CUT_MEDIAN_FLUX and med_flux > ABS_MIN:
        if MIN_W == 0:
            MIN_W=min_wave
        if MAX_W == 0:
            MAX_W=max_wave
        ################
        # REDSHIFT DETERMINATION   
        my_plot=plot
        
        fit_redshift=0
        redshift_abs=redshift
        delta_chi=10
        NITER=0
        niter_tmp_max=10
        chi_sq_min_now=1e12
        min_chi_sq_limit=min_chi_sq
        n_mc=20
        plot=my_plot
        #pdl_rat_master=np.ones(n_unc+1)
        [min_chi_sq,pdl_age_mod,pdl_met_mod,pdl_ml,pdl_Av,coeffs,coeffs_N,coeffs_NM,pdl_model_spec_min,pdl_res]=ssp1.fit_ssp_lin_no_zero(redshift,sigma,Av,crval_kin,cdelt_kin,crpix_kin,nf_kin,n_c_kin,pdl_flux_c_ini_kin,hdr2,wave_unc,masked_Av,e_flux_unc,flux_unc,n_mc,chi_sq_min_now,min_chi_sq,plot,sigma_inst)
        #
        # We substract the continuum!
        #
        #sys.exit()
        n_mc=40
        pdl_mod_SSP=pdl_model_spec_min
        pdl_res_SSP=pdl_res
        i0_now=int(0.4*nx)
        i1_now=int(0.6*nx)
        #print pdl_res[i0_now:i1_now],rank
        #print pdl_model_spec_min[i0_now:i1_now],rank
        #sys.exit()
        stats_res=np.std(pdl_res[i0_now:i1_now])+np.mean(pdl_res[i0_now:i1_now])
        stats_mod=np.mean(pdl_model_spec_min[i0_now:i1_now])
        SN=0
        if stats_res > 0:
            SN=abs(stats_mod/stats_res)
        print "Signal-to-Noise = "+str(SN)+" ("+str(stats_mod)+"/"+str(stats_res)+") at processor "+str(rank)
        old=1       
        if old == 1 and SN > 10:
            pdl_model_spec_min[np.where(pdl_model_spec_min == 0)[0]]=1.
            pdl_rat=pdl_res/pdl_model_spec_min+1
            rat=pdl_rat
            sigma_mean=np.sqrt(sigma_inst**2.0+(5000.0*(sigma/300000.0))**2.0)
            med_rat=my.median_filter(int(5*2.354*sigma_mean),rat)
            pdl_med_rat=np.array(med_rat)
            n_unc_1=n_unc-1
            pdl_wave_unc=wave_unc[0]+(wave_unc[1]-wave_unc[0])*np.arange(0,n_unc_1)
            med_rat=my.median_filter(int(7*2.354*sigma_mean),rat)
            med_sigma=int(1.5*sigma)
            if med_sigma < 3:
                med_sigma=3
            med_rat_box=my.median_box(med_sigma,med_rat)
            med_wave_box=my.median_box(med_sigma,wave_unc)
            pdl_wave_unc=np.array(wave_unc)
            pdl_wave_box=np.array(med_wave_box)
            pdl_med_rat_box=np.array(med_rat_box)
            DD1=pdl_wave_unc.shape
            DD2=pdl_wave_box.shape
            DD3=pdl_med_rat_box.shape
            if DD2[0] > DD3[0]:
                ND=DD3[0]-1
                pdl_wave_box=pdl_wave_box[0:ND]
            else:
                if DD2[0] < DD3[0]:
                    ND=DD2[0]-1
                    pdl_med_rat_box=pdl_med_rat_box[0:ND] 
            y_rat = interp1d(pdl_wave_box, pdl_med_rat_box,bounds_error=False,fill_value=0.)(pdl_wave_unc)
            if plot > 0:
                out_ps_now="junk2"
                title="ratio"
                ssp.plot_results_min_max(2,wave_unc,[flux_unc,pdl_model_spec_min,pdl_res,pdl_rat,y_rat],out_ps_now,title,-0.5,2.0)
            i0_now=int(0.4*nx)
            i1_now=int(0.6*nx)
            stats_rat=my.stats(y_rat[i0_now:i1_now])
            if stats_rat[0] > 0.8 and stats_rat[0] < 1.2 and stats_rat[1] > 0.02:
                for i in range(0, n_unc):
                    val=y_rat[i]
                    if val > (stats_rat[3]-stats_rat[1]) and val < (stats_rat[3]+stats_rat[1]) and val != 0:
                        flux_unc[i]=flux_unc[i]/val
                        flux_unc_org[i]=flux_unc_org[i]/val
        coeffs_cat=np.zeros([nf+1,n_mc])
        NITER=0
        while MIN_CHISQ > MIN_DELTA_CHISQ and NITER < MAX_NITER:
            if NITER == 1:
                MIN_CHISQ=1e12
            #print "("+str(MIN_CHISQ)+">"+str(MIN_DELTA_CHISQ)+")&&("+str(NITER)+"<"+str(MAX_NITER)+") "+str(med_flux)+">"+str(CUT_MEDIAN_FLUX)+" "+str(jy)+"/"+str(NY)
            if med_flux > CUT_MEDIAN_FLUX:
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
                    #f=open(out_elines,"a")
                    #f.write("#ID "+str(jy)+"\n")
                    #f.close
                    e_line_jy.extend(["#ID "+str(jy)+"\n"])
                    for ist in range(0,ns):
                        SYS_VEL_IN=input_redshift*vel_light
                        SYS_VEL_MIN_IN=input_min_redshift*vel_light
                        SYS_VEL_MAX_IN=input_max_redshift*vel_light
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
                        if SYS_VEL < SYS_VEL_MIN_IN or SYS_VEL > SYS_VEL_MAX_IN:
                            SYS_VEL = SYS_VEL_IN
                        if SYS_VEL_MAX > SYS_VEL_MAX_IN:
                            SYS_VEL_MAX = SYS_VEL_MAX_IN
                        if SYS_VEL_MIN < SYS_VEL_MIN_IN:
                            SYS_VEL_MIN=SYS_VEL_MIN_IN
                        start_w_e=start_w_E[ist]
                        end_w_e=end_w_E[ist]
                        mask_e=mask_E[ist]
                        config_e=config_E[ist]
                        npoly_e=npoly_E[ist]
                        mask_poly_e=mask_poly_E[ist]
                        nmin_e=nmin_E[ist]
                        nmax_e=nmax_E[ist]   
                        #print "CONF="+config_e
                        wave_elines=[]
                        flux_elines=[]
                        flux_e_elines=[]
                        masked_elines=[]
                        n_e=0
                        for i in range(0, n_unc):
                            if wave_unc[i] > start_w_e and wave_unc[i] < end_w_e:
                                wave_elines.extend([wave_unc[i]])
                                flux_elines.extend([flux_unc_org[i]-pdl_mod_SSP[i]])
                                flux_e_elines.extend([e_flux_unc[i]])
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
                            d_time=ssp.get_seconds()-time
                            chi_single=0
                            #print np.sum(pdl_flux_e_elines)
                            #[chi_sq_now,pdl_a,pdl_model_tmp,pdl_model_cont_tmp]=ssp.fit_elines_grad_rnd_new(pdl_wave_elines,pdl_flux_elines,pdl_flux_e_elines,n_mod,chi_goal,d_chi_goal,typef,a_out_tmp,ia,a0,a1,link,n_mc_e,pdl_masked_elines,deft,scale_ini)#,max_time)
                            [chi_sq_now,pdl_a,pdl_model_tmp,pdl_model_cont_tmp]=ssp.fit_elines_grad_rnd_new(pdl_wave_elines,pdl_flux_elines,pdl_flux_e_elines,n_mod,chi_goal,d_chi_goal,typef,a_out_tmp,ia,a0,a1,link,n_mc_e,pdl_masked_elines,deft,scale_ini)#,max_time)
                            a_out_now=ssp.copy_a_pdl(n_mod,pdl_a)
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
                        if len(pdl_model.shape) == 1:
                            arp1=pdl_model
                        else:
                            arp1=pdl_model[0,:]
                        ssp.plot_results(plot,pdl_wave_elines,[pdl_flux_elines,arp1,(pdl_flux_elines-arp1)],out_ps_now,title)
                        #print "----------------------------------------";
                        a_final=ssp.mean_a_results_last(n_mod,nnk,a_results,ia)
                        pdl_res_now=pdl_flux_elines-pdl_model
                        stats_back=my.stats(pdl_res_now)
                        a_final=ssp.add_back_noise(n_mod,a_final,typef,chi_sq_now,stats_back[1])                        
                        #ssp.print_a_final(n_mod,a_final,typef,chi_sq_now)#Only for test
                        out_fit_spectra=out_elines
                        temp_f=ssp.print_a_final_file_add_mpi(n_mod,a_final,typef,chi_sq_now,out_fit_spectra,jy)
                        e_line_jy.extend(temp_f)
                        SYS_VEL_elines=a_final[0][0][3]
                        red_elines=SYS_VEL_elines/vel_light  
                        [n_mod_fixed,junk_a_fixed,junk_a_type_fixed]=ssp.add_a_results_elines(n_mod,a_final,typef,n_mod_fixed,a_fixed,a_type_fixed)
                        a_fixed=junk_a_fixed
                        a_type_fixed=junk_a_type_fixed
                        nmin_e=int(0.1*n_unc)
                        nmax_e=int(0.9*n_unc)    
                        ###############################
                        # Low order polynomical!
                        out_fit_now=out_fit+"."+str(start_w_e)+"_"+str(end_w_e)+".pdf"
                        sigma_mean=np.sqrt(sigma_inst**2.0+(5000.0*(sigma/300000.0))**2.0)
                        box=int(sigma_mean*6)
                        #print "DONE FIT ELINES CONFIG "+str(ist)
                        #
                        # We create a FIXED model of the emission lines
                        #
                    pdl_model_elines=np.zeros(n_unc)
                    pdl_model_cont=np.zeros(n_unc)
                    pdl_wave_elines=np.array(wave_unc)
                    i_loops_ssp=0
                    for i in range(0, n_mod_fixed):
                        pdl_tmp=ssp.create_single_model(pdl_wave_elines,i,a_type_fixed,a_fixed)
                        NN=n_unc
                        pdl_model_elines=pdl_model_elines+pdl_tmp[0,:]
                    #
                    # We remove the gas before a new iteration
                    #
                    for i in range(0, n_unc):
                        flux_unc[i]=flux_unc_org[i]-pdl_model_elines[i]
                    pdl_mod_JOINT=pdl_mod_SSP+pdl_model_elines
                    pdl_res_JOINT=pdl_res_SSP-pdl_model_elines
                    pdl_no_gas=np.array(flux_unc)
                    ##############################################################
                    # We rescale!
                    ##############################################################
                    y_rat=np.ones(nx+1)
                    if SN > 10:
                        pdl_mod_JOINT[np.where(pdl_mod_JOINT == 0)[0]]=1.
                        pdl_rat=pdl_res_JOINT/pdl_mod_JOINT+1
                        rat=pdl_rat
                        n_unc_1=n_unc-1
                        pdl_wave_unc=wave_unc[0]+(wave_unc[1]-wave_unc[0])*np.arange(0, n_unc_1)
                        sigma_mean=np.sqrt(sigma_inst**2.0+(5000.0*(sigma/300000.0))**2.0)
                        med_rat=my.median_filter(int(5*2.354*sigma_mean),rat)
                        med_sigma=int(1.5*sigma_mean)
                        if med_sigma < 3:
                            med_sigma=3
                        med_rat_box=my.median_box(med_sigma,med_rat)
                        med_wave_box=my.median_box(med_sigma,wave_unc)
                        pdl_wave_unc=np.array(wave_unc)
                        pdl_wave_box=np.array(med_wave_box)
                        pdl_med_rat_box=np.array(med_rat_box)
                        DD1=pdl_wave_unc.shape
                        DD2=pdl_wave_box.shape
                        DD3=pdl_med_rat_box.shape
                        if DD2[0] > DD3[0]:
                            ND=DD3[0]-1
                            pdl_wave_box=pdl_wave_box[0:ND]
                        else:
                            if DD2[0] < DD3[0]:
                                ND=DD2[0]-1
                                pdl_med_rat_box=pdl_med_rat_box[0:ND]
                        y_rat = interp1d(pdl_wave_box, pdl_med_rat_box,bounds_error=False,fill_value=0.)(pdl_wave_unc) 
                        i0_now=int(0.4*nx)
                        i1_now=int(0.6*nx)
                        stats_rat0=np.mean(y_rat[i0_now:i1_now])
                        stats_rat1=np.std(y_rat[i0_now:i1_now])+stats_rat0
                        stats_rat2=np.median(y_rat[i0_now:i1_now])
                        stats_rat3=np.amin(y_rat[i0_now:i1_now])        
                        if plot > 1:
                            out_ps_now="junk3"
                            title="ratio = "+str(stats_rat0)+", rms="+str(stats_rat1)
                            print title
                            ssp.plot_results_min_max(2,wave_unc,[flux_unc,pdl_model_spec_min,pdl_res,pdl_rat,y_rat],out_ps_now,title,0.0,2.0)
                        if stats_rat0 > 0.8 and stats_rat0 > 1.2 and stats_rat1 > 0.02:
                            if jy == 0:
                                # Continuum shape correction on/off
                                pdl_rat_master=y_rat
                                pdl_rat_master[np.where(pdl_rat_master == 0)[0]]=1.
                            else:
                                y_rat=pdl_rat_master
                        else:
                            y_rat=pdl_rat_master
                    else:
                        i0_now=int(0.4*nx)
                        i1_now=int(0.6*nx)
                        stats_rat0=np.mean(y_rat[i0_now:i1_now])
                        stats_rat1=np.std(y_rat[i0_now:i1_now])+stats_rat0
                        stats_rat2=np.median(y_rat[i0_now:i1_now])
                        stats_rat3=np.amin(y_rat[i0_now:i1_now])
                    for i in range(0, n_unc):
                        val=y_rat[i]
                        if val > (stats_rat3-stats_rat1) and val < (stats_rat3+stats_rat1):
                            flux_unc[i]=flux_unc[i]/val
                            flux_unc_org[i]=flux_unc_org[i]/val
                            ##############################################################
                            # End re-scale
                            ##############################################################
                    #pdl_res_SSP_now=flux_unc_org-pdl_mod_SSP
                    #statst=my.stats(pdl_res_SSP_now)                            
                    #print "Deriving SFH.... "
                    #########################################################
                    # START Fitting the SFH
                    ######################################################### 
                    pdl_e_flux_now=e_flux_unc+np.sqrt(abs(pdl_model_elines))
                    e_flux_now=pdl_e_flux_now
                    n_mc=10
                    #print "Start full SFH"
                    sigma_MaNGA=np.sqrt(abs(sigma**2.-2.**2.))                     
                    [min_chi_sq,pdl_age_mod,pdl_met_mod,pdl_ml,pdl_Av,coeffs,coeffs_N,coeffs_NM,pdl_mod_SSP,pdl_res_SSP,coeffs_N_input,e_coeffs_N_input]=ssp1.fit_ssp_lin_MC(redshift,sigma,Av,crval1,cdelt1,crpix1,nf,n_c,pdl_flux_c_ini,hdr,wave_unc,masked,e_flux_unc,flux_unc,n_mc,chi_sq_min_now,MIN_CHISQ,plot,sigma_inst)
                    sigma_mean=np.sqrt(sigma_inst**2.0+(5000.0*(sigma/300000.0))**2.0)
                    smooth_ratiot=ssp.smooth_ratio(flux_unc,pdl_mod_SSP,int(sigma_mean))
                    pdl_mod_SSP_no_cor=np.copy(pdl_mod_SSP)
                    pdl_mod_SSP=pdl_mod_SSP*smooth_ratiot
                    #print "End full SFH"
                    i_loops_ssp=i_loops_ssp+1
                    pdl_mod_JOINT=pdl_mod_SSP+pdl_model_elines
                    pdl_res_JOINT=pdl_res_SSP-pdl_model_elines
                    print "SSP From processor "+str(rank)
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
                    a_coeffs_N=coeffs_N[:,0]  
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
                    ML=0
                    sum_ML=0
                    sum_I=0
                    sum_II=0
                    for k in range(0, nf):
                        if a_coeffs[k] > 0:
                            a_e_coeffs_N[k]=a_e_coeffs[k]*(a_coeffs_N[k]/a_coeffs[k])
                        else:
                            a_e_coeffs_N[k]=0
                        l_age_min=l_age_min+a_min_coeffs[k]*np.log10(age_mod[k])
                        l_met_min=l_met_min+a_min_coeffs[k]*np.log10(met_mod[k])
                        l_Av_min=l_Av_min+a_min_coeffs[k]*np.log10(Av_mod[k])
                        l_age_min_mass=l_age_min_mass+(ml[k])*a_coeffs_N[k]*np.log10(age_mod[k])
                        l_met_min_mass=l_met_min_mass+(ml[k])*a_coeffs_N[k]*np.log10(met_mod[k])
                        l_Av_min_mass=l_Av_min_mass+(ml[k])*a_coeffs_N[k]*np.log10(Av_mod[k])
                        e_l_age_min=e_l_age_min+a_e_coeffs[k]*np.log10(age_mod[k])
                        e_l_met_min=e_l_met_min+a_e_coeffs[k]*np.log10(met_mod[k])
                        e_l_Av_min=e_l_Av_min+a_e_coeffs[k]*np.log10(Av_mod[k])                
                        e_l_age_min_mass=e_l_age_min_mass+(ml[k])*a_e_coeffs_N[k]*np.log10(age_mod[k])
                        e_l_met_min_mass=e_l_met_min_mass+(ml[k])*a_e_coeffs_N[k]*np.log10(met_mod[k])
                        e_l_Av_min_mass=e_l_Av_min_mass+(ml[k])*a_e_coeffs_N[k]*np.log10(Av_mod[k])
                        ML=ML+ml[k]*a_coeffs_N[k]
                        sum_ML=sum_ML+ml[k]*a_coeffs_N[k]
                        sum_I=sum_I+a_coeffs_N[k]
                        sum_II=sum_II+a_min_coeffs[k]
                        guess=coeffs_N_input[k]
                        #f=open(out_coeffs_file,"a")
                        #f.write(("%2d" % k)+" "+("%7.4f" % age_mod[k])+" "+("%7.4f" % met_mod[k])+" "+("%7.4f" % a_coeffs_N[k])+" "+("%7.4f" % a_min_coeffs[k])+" "+("%4.4f" % np.log10(ml[k]))+" "+("%4.2f" % Av_mod[k])+" "+("%7.4f" % a_coeffs[k])+" "+("%7.4f" % a_e_coeffs[k])+" "+("%3d" % jy)+" "+("%7.4f" % guess)+"\n")
                        #f.close
                        out_coeffs_f.extend([("%2d" % k)+" "+("%7.4f" % age_mod[k])+" "+("%7.4f" % met_mod[k])+" "+("%7.4f" % a_coeffs_N[k])+" "+("%7.4f" % a_min_coeffs[k])+" "+("%4.4f" % np.log10(ml[k]))+" "+("%4.2f" % Av_mod[k])+" "+("%7.4f" % a_coeffs[k])+" "+("%7.4f" % a_e_coeffs[k])+" "+("%3d" % jy)+" "+("%7.4f" % guess)+"\n"])
                        if a_coeffs[k] > 1e-5:
                            print ("%2d" % k)+" "+("%7.4f" % age_mod[k])+" "+("%7.4f" % met_mod[k])+" "+("%7.4f" % a_coeffs_N[k])+" "+("%7.4f" % a_min_coeffs[k])+" "+("%4.4f" % np.log10(ml[k]))+" "+("%4.2f" % Av_mod[k])+" "+("%7.4f" % a_coeffs[k])+" "+("%7.4f" % a_e_coeffs[k])+" "+("%7.4f" % guess)
                    print "------------------------------------------------------------------------------"
                    #print "SUM_I = "+str(sum_I)
                    #print "SUM_II = "+str(sum_II) 
                    age_min=10**(l_age_min)
                    met_min=10**(l_met_min)
                    Av_min=10**(l_Av_min)
                    if ML > 0:
                        age_min_mass=10**(l_age_min_mass/sum_ML)
                        met_min_mass=10**(l_met_min_mass/sum_ML)
                        Av_min_mass=10**(l_Av_min_mass/sum_ML)
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
                #####################################################################
                # END fitting SFH
                #####################################################################
                # Copy output!
                pdl_output[jy,:,0]=flux_unc_org
                pdl_output[jy,:,1]=pdl_mod_SSP
                pdl_output[jy,:,2]=pdl_mod_JOINT
                pdl_res_SSP=np.array(flux_unc_org)-pdl_mod_SSP
                pdl_res_SSP_no_cor=np.array(flux_unc_input)-pdl_mod_SSP_no_cor
                pdl_output[jy,:,3]=pdl_res_SSP_no_cor
                pdl_tmp=np.array(flux_unc_org)
                pdl_res_JOINT=pdl_tmp-pdl_mod_JOINT
                pdl_output[jy,:,4]=pdl_res_JOINT
                # NO GAS    
                pdl_output[jy,:,5]=np.array(flux_unc_org)-(pdl_res_SSP-pdl_res_JOINT)
                pdl_output_now[:,0]=np.array(flux_unc_org)
                pdl_output_now[:,1]=pdl_mod_SSP
                pdl_output_now[:,2]=pdl_mod_JOINT
                pdl_output_now[:,3]=pdl_res_SSP
                pdl_output_now[:,4]=pdl_res_JOINT
                pdl_output_now[:,5]=np.array(flux_unc_org)-(pdl_res_SSP-pdl_res_JOINT)
                title="X="+str(chi_sq_now)+" T="+str(age_min)+" ("+str(age_min_mass)+") Z="+str(met_min)+" ("+str(met_min_mass)+") Av="+str(Av_min)+" z="+str(redshift)+" sigma="+str(sigma)
            #print "I.Iter = "+str(NITER)+" DONE"
            NITER=NITER+1    
        print "--------------------------------------------------------------" 
        pdl_masked=np.array(masked)
        pdl_chi_now=((pdl_masked*pdl_res_JOINT)**2)/((np.array(e_flux_unc))**4)
        pdl_chi_now[np.isnan(pdl_chi_now)]=0
        chi_joint=np.sum(pdl_chi_now)
        chi_joint=(chi_joint/(n_unc-n_mod_fixed-nf-1))**0.5
        rms=np.std(np.array(pdl_masked*pdl_res_JOINT))
        FLUX=0
    else:
        [chi_joint,age_min,e_age_min,met_min,e_met_min,Av_min,e_Av_min,out_redshift,e_redshift,out_sigma,e_sigma,FLUX,redshift_abs,med_flux,rms,age_min_mass,met_min_mass,SYS_VEL]=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
        #f=open(out_coeffs_file,"a")
        for k in range(0, nf):
            #f.write(("%2d" % k)+" "+("%7.4f" % age_mod[k])+" "+("%7.4f" % met_mod[k])+" "+("%7.4f" % a_coeffs_N[k])+" "+("%7.4f" % a_min_coeffs[k])+" "+("%4.4f" % np.log10(ml[k]))+" "+("%4.2f" % Av_mod[k])+" "+("%7.4f" % a_coeffs[k])+" "+("%7.4f" % a_e_coeffs[k])+" "+("%3d" % jy)+" "+("%7.4f" % guess)+"\n")
            out_coeffs_f.extend([("%2d" % k)+" "+("%7.4f" % age_mod[k])+" "+("%7.4f" % met_mod[k])+" "+("%7.4f" % a_coeffs_N[k])+" "+("%7.4f" % a_min_coeffs[k])+" "+("%4.4f" % np.log10(ml[k]))+" "+("%4.2f" % Av_mod[k])+" "+("%7.4f" % a_coeffs[k])+" "+("%7.4f" % a_e_coeffs[k])+" "+("%3d" % jy)+" "+("%7.4f" % guess)+"\n"])
        #f.close        
    title="X="+str(chi_sq_now)+" T="+str(age_min)+" ("+str(age_min_mass)+") Z="+str(met_min)+" ("+str(met_min_mass)+") Av="+str(Av_min)+" z="+str(redshift)+" sigma="+str(sigma)
    if pdl_output_now[:,0:4].shape[0] == len(wave_unc):
        pdl_output_temp=np.transpose(pdl_output_now[:,0:4])
    else:
        pdl_output_temp=pdl_output_now[:,0:4]    
    ssp.plot_results_min_max(plot,wave_unc,pdl_output_temp,out_ps,title,y_min,y_max)
    ZH=np.log10(met_min/0.02)
    lML=np.log10(ML)
    print "MSP CHISQ="+str(chi_joint)+" AGE="+str(age_min)+"+-"+str(e_age_min)+" MET="+str(met_min)+"+-"+str(e_met_min)+" AV="+str(Av_min)+"+-"+str(e_Av_min)+" REDSHIFT="+str(redshift)+"+-"+str(e_redshift)+" SIGMA_DISP="+str(sigma)+"+-"+str(e_sigma)+" RMS="+str(rms)+" MED_FLUX="+str(med_flux)+" AGE_mass="+str(age_min_mass)+"+-"+str(e_age_min_mass)+" MET_mass="+str(met_min_mass)+"+-"+str(e_met_min_mass)+" ZH_lum="+str(ZH)+" log_M/L="+str(lML)
    j1=int(0.4*n_unc)
    j2=int(0.6*n_unc)
    wave_norm=(wave_unc[j1]+wave_unc[j2])/2.
    if chi_joint == 0:
        chi_joint=1
    #f=open(outfile, "a")
    #f.write(str(chi_joint)+","+str(age_min)+","+str(e_age_min)+","+str(met_min)+","+str(e_met_min)+","+str(Av_min)+","+str(e_Av_min)+","+str(redshift)+","+str(e_redshift)+","+str(sigma)+","+str(e_sigma)+","+str(FLUX)+","+str(redshift_abs)+","+str(med_flux)+","+str(rms)+","+str(age_min_mass)+","+str(e_age_min_mass)+","+str(met_min_mass)+","+str(e_met_min_mass)+","+str(SYS_VEL)+","+str(lML)+","+str(jy)+"\n")
    #f.close
    openfile_f.extend([str(chi_joint)+","+str(age_min)+","+str(e_age_min)+","+str(met_min)+","+str(e_met_min)+","+str(Av_min)+","+str(e_Av_min)+","+str(redshift)+","+str(e_redshift)+","+str(sigma)+","+str(e_sigma)+","+str(FLUX)+","+str(redshift_abs)+","+str(med_flux)+","+str(rms)+","+str(age_min_mass)+","+str(e_age_min_mass)+","+str(met_min_mass)+","+str(e_met_min_mass)+","+str(SYS_VEL)+","+str(lML)+","+str(jy)+"\n"])
    print str(jy)+"/"+str(NY)+" DONE from processor "+str(rank)
    print "--------------------------------------------------------------" 
#
# End of the RSS loop!
#
# Write output file
if rank != 0:
    comm.send(openfile_f, dest=0, tag=13)
    comm.send(out_coeffs_f, dest=0,tag=12)
    comm.send(e_line_jy, dest=0,tag=11)
    comm.send(pdl_output, dest=0,tag=10)
if rank == 0:
    for i in range(1, nprocs):
        open_ft=comm.recv(source=i, tag=13)
        openfile_f.extend(open_ft)
        out_ft=comm.recv(source=i, tag=12)
        out_coeffs_f.extend(out_ft)
        elin_ft=comm.recv(source=i, tag=11)
        e_line_jy.extend(elin_ft)
        pdl_outputt=comm.recv(source=i, tag=10)
        i1=ind1[i]
        i2=ind2[i]
        pdl_output[i1:i2,:]=pdl_outputt[i1:i2,:]
    f=open(outfile, "a")
    for i in range(0, len(openfile_f)):
        f.write(openfile_f[i])
    f.close
    f=open(out_coeffs_file, "a")
    for i in range(0, len(out_coeffs_f)):
        f.write(out_coeffs_f[i])
    f.close
    f=open(out_elines, "a")
    for i in range(0, len(e_line_jy)):
        f.write(e_line_jy[i])
    f.close
    h=pyf.PrimaryHDU().header
    h["NAXIS"]=3
    h["NAXIS1"]=nx 
    h["NAXIS2"]=NY
    h["NAXIS3"]=6
    h["COMMENT"]="output FITS file"
    h['NAME0']="org_spec"
    h['NAME1'] ="model_spec"
    h['NAME2'] ="mod_joint_spec"
    h['NAME3']="gas_spec"
    h['NAME4']="res_joint_spec"
    h['NAME5']="no_gas_spec"
    h["CRVAL1"]=crval
    h["CDELT1"]=cdelt
    h["CRPIX1"]=crpix
    h['FILENAME']=out_fit
    if pt.exists(out_fit) == False:
        wfits(out_fit,pdl_output,h)
    else:
        sycall("rm "+out_fit,1)
        wfits(out_fit,pdl_output,h)
    ################################
    sec_end=ssp.print_time()
    sec_total=sec_end-sec_ini
    print "# SECONDS = "+str(sec_total)