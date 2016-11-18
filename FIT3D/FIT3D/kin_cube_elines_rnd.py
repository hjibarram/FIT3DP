#!/usr/bin/python
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
from pyfits import writeto as wfit
from scipy.interpolate.interpolate import interp1d
import my as my
import os.path as ptt
import matplotlib
matplotlib.use('Agg')

def sycall(comand):
    from subprocess import call
    line=comand.split(" ")
    fcomand=[]
    fcomand.extend(line)
    call(fcomand)

def wfits(name, data, hdr):
    if ptt.exists(name) == False:
        wfit(name,data,hdr)
    else:
        sycall("rm "+name)
        wfit(name,data,hdr)
        
def mean_a_results(n_mod,nk,a_print,I,J,type,pdl_wave):
    a_copy=np.zeros([2,n_mod,9])
    for ii in range(0, n_mod):
        for j in range(0, 9):
            a_tmp=[]
            kk=0
            val=a_print[0][ii][j]
            for k in range(1, nk):
                a_tmp.extend([a_print[k][ii][j]])
                kk=kk+1
            if kk > 0:
                if kk > 3:
                    val=np.median(a_tmp)
                else:
                    val=np.mean(a_tmp)
                e_val=np.std(a_tmp)
            else:
                val=0
                e_val=0
            a_copy[0][ii][j]=val
            a_copy[1][ii][j]=e_val*ia[ii][j]
        if type[ii] == "poly1d\n":
            sum=0;
            sum_e=0;
            a_wave=pdl_wave
            for jj in range(0, len(a_wave)):
                for j in range(0, 9):
                    sum=sum+a_copy[0][ii][j]*(a_wave[jj])**j
                    sum_e=sum_e+a_copy[1][ii][j]*(a_wave[jj])**j
            pdl_f[ii,I,J]=sum
            pdl_ef[ii,I,J]=sum_e
        else:        
            pdl_f[ii,I,J]=a_copy[0][ii][1]
            pdl_ef[ii,I,J]=a_copy[1][ii][1]
            pdl_dw[ii,I,J]=2.354*a_copy[0][ii][2]
            pdl_edw[ii,I,J]=2.354*a_copy[1][ii][2]
            pdl_dw_s[ii,I,J]=a_copy[0][ii][2]
            pdl_edw_s[ii,I,J]=a_copy[1][ii][2]
            pdl_v[ii,I,J]=a_copy[0][ii][3]
            pdl_ev[ii,I,J]=a_copy[1][ii][3]
    return a_copy

        
sys.argv=filter(None,sys.argv)
if len(sys.argv) < 7:
    print "USE: kin_cube_elines_rnd.py CUBE.fits CONFIG MASK_LIST W_MIN W_MAX OUTFILE [N_MC=50] [N_LOOPS=5] [PLOT] [SCALE_INI=0.15] [PREFIX] [MEMO=0/1] [input_vel_map.fits,input_vel_mask_now,FIXED=0,1] [input_disp_map.fits,FIXED=0,1 | FWHM=2.354*sigma]"
    sys.exit(0)
sec_ini=ssp.print_time()
input_cube=sys.argv[1]
config=sys.argv[2]
mask_list=sys.argv[3]
w_min=float_(sys.argv[4])
w_max=float_(sys.argv[5])
out_fit_spectra=sys.argv[6]
prefix="map_"+out_fit_spectra

n_loops=3
n_mc=50
deft=0

if len(sys.argv) == 8:
    n_mc=int(sys.argv[7])

if len(sys.argv) == 9:
    n_mc=int(sys.argv[7])
    n_loops=int(sys.argv[8])
plot=0

if len(sys.argv) == 10:
    n_mc=int(sys.argv[7])
    n_loops=int(sys.argv[8])
    plot=int(sys.argv[9])

scale_ini=0.15
if len(sys.argv) == 11:
    n_mc=int(sys.argv[7])
    n_loops=int(sys.argv[8])
    plot=int(sys.argv[9])
    scale_ini=float_(sys.argv[10])

if len(sys.argv) == 12:
    n_mc=int(sys.argv[7])
    n_loops=int(sys.argv[8])
    plot=int(sys.argv[9])
    scale_ini=float_(sys.argv[10])
    prefix=sys.argv[11]

memo=0
if len(sys.argv) == 13:
    n_mc=int(sys.argv[7])
    n_loops=int(sys.argv[8])
    plot=int(sys.argv[9])
    scale_ini=float_(sys.argv[10])
    prefix=sys.argv[11]
    memo=float_(sys.argv[12])

guided=0
guided_type=2
if len(sys.argv) == 14:
    n_mc=int(sys.argv[7])
    n_loops=int(sys.argv[8])
    plot=int(sys.argv[9])
    scale_ini=float_(sys.argv[10])
    prefix=sys.argv[11]
    memo=float_(sys.argv[12])
    JUNK=sys.argv[13]
    guided=1
    [guided_map,guided_mask,guided_type]=JUNK.split(",")
    gmap=gdata(guided_map)
    gmask=gdata(guided_mask)
    guided_type=float_(guided_type)    

guided_disp=0
guided_disp_type=2
if len(sys.argv) == 15:
    n_mc=int(sys.argv[7])
    n_loops=int(sys.argv[8])
    plot=int(sys.argv[9])
    scale_ini=float_(sys.argv[10])
    prefix=sys.argv[11]
    memo=float_(sys.argv[12])
    JUNK=sys.argv[13]
    guided=1
    [guided_map,guided_mask,guided_type]=JUNK.split(",")
    gmap=gdata(guided_map)
    gmask=gdata(guided_mask)    
    JUNK=sys.argv[14]
    guided_disp=1
    [guided_disp_map,guided_disp_type]=JUNK.split(",")
    gmap_disp=gdata(guided_disp_map)
    gmap_disp=gmap_disp/2.354;
    stats_disp=my.stats(gmap_disp)
    guided_type=float_(guided_type)
    guided_disp_type=float_(guided_disp_type)    

##################################################
# READING INPUTS
##################################################

#
# Reading MASK LIST
#
nmask=0
start_mask=[]
end_mask=[]
if mask_list == "none":
    nmask=0
else:
    f=open(mask_list, 'r')
    for line in f:
        data=line.split(" ")
        data=filter(None,data)
        start_mask.extend([float_(data[0])])
        end_mask.extend([float_(data[1])])
        nmask=nmask+1
    f.close()

print str(nmask)+" regions to mask"
#
# READING FILE!
#
[pdl_input_cube, h]=gdata(input_cube, 0, header=True)
[nz,nx,ny]=pdl_input_cube.shape
crpix=h["CRPIX3"]
crval=h["CRVAL3"]
cdelt=h["CDELT3"]

if w_min < crval:
    w_min=crval

if w_max >= (crval+cdelt*nz):
    w_max=(crval+cdelt*(nz-1))

pdl_wave_org=crval+cdelt*np.arange(0,nz)
i0=int((w_min-crval)/cdelt)
i1=int((w_max-crval)/cdelt)
pdl_flux_cube=pdl_input_cube[i0:i1,:,:]
pdl_wave=pdl_wave_org[i0:i1]
nm=i1-i0
pdl_masked=np.ones(nz)
org=np.zeros([nm,nx,ny])
mod=np.zeros([nm,nx,ny])
res=np.zeros([nm,nx,ny])

for i in range(0, nz):
    w_test=pdl_wave_org[i]
    pdl_masked[i]=1
    for j in range(0, nmask):
        if w_test > start_mask[j] and w_test < end_mask[j]:
            pdl_masked[i]=0
pdl_masked_org=pdl_masked
pdl_masked=pdl_masked_org[i0:i1]

#
# Reading CONFIG!
#
typef=[]
f=open(config,"r")
line=f.readline()
data=line.split(" ")
data=filter(None,data)
[junk,n_mod,chi_goal,d_chi_goal]=data
n_mod=int(n_mod)
a=np.zeros([n_mod,9])
ia=np.zeros([n_mod,9])
ea=np.zeros([n_mod,9])
a0=np.zeros([n_mod,9])
a1=np.zeros([n_mod,9])
link=np.zeros([n_mod,9])
a_IN=np.zeros([n_mod,9])
ia_IN=np.zeros([n_mod,9])
a0_IN=np.zeros([n_mod,9])
a1_IN=np.zeros([n_mod,9])
link_IN=np.zeros([n_mod,9])
for ii in range(0, n_mod):
    line=f.readline()
    data=line.split(" ")
    if data[0] == "eline\n" or data[0] == "poly1d\n":
        cnf=data[0]
        typef.extend([cnf])
    for j in range(0, 9):
        line=f.readline()
        data=line.split(" ")
        data=filter(None,data)
        data=[g for g in data if g != "\t"]
        data0=[]
        for ta in data:
            temp=ta.split('\t')
            for tb in temp:
                data0.extend([tb])
        data=data0
        data=filter(None,data)
        a[ii][j]=float_(data[0])
        ia[ii][j]=float_(data[1])
        ea[ii][j]=0
        a0[ii][j]=float_(data[2])
        a1[ii][j]=float_(data[3])
        link[ii][j]=float_(data[4])
        a_IN[ii][j]=float_(data[0])
        ia_IN[ii][j]=float_(data[1])
        a0_IN[ii][j]=float_(data[2])
        a1_IN[ii][j]=float_(data[3])
        link_IN[ii][j]=float_(data[4])
f.close()
print str(n_mod)+" models to fit ("+str(n_mod)+","+str(nx)+","+str(ny)+")"

pdl_f=np.zeros([n_mod,nx,ny])
pdl_ef=np.zeros([n_mod,nx,ny])
pdl_v=np.zeros([n_mod,nx,ny])
pdl_ev=np.zeros([n_mod,nx,ny])
pdl_dw=np.zeros([n_mod,nx,ny])
pdl_edw=np.zeros([n_mod,nx,ny])
pdl_dw_s=np.zeros([n_mod,nx,ny])
pdl_edw_s=np.zeros([n_mod,nx,ny])

for ix in range(0, nx):
    for jy in range(0, ny):
        pdl_flux=pdl_flux_cube[:,ix,jy]
        pdl_flux[np.isnan(pdl_flux)]=0 
        statst=my.stats(pdl_flux)
        val_median=abs(statst[2])    
        val_sigma=abs(statst[1])
        pdl_e_flux=(val_median+val_sigma**2.0)*np.ones(nm)
        if val_median != 0:
            pdl_e_flux[np.where(pdl_e_flux == 0)[0]]=val_median
        else:
            pdl_e_flux[np.where(pdl_e_flux == 0)[0]]=0.001
        pdl_model=np.zeros(nm)
        pdl_model_cont=np.zeros(nm)
        chi_sq_now=5
        chi_sq_start=chi_sq_now
        a_results=np.zeros([1, n_mod, 9])
        nnk=0
        i_ter=0
        a_out_now=ssp.copy_a(n_mod,a)
        ia_now=ssp.copy_a(n_mod,ia)
        a0_now=ssp.copy_a(n_mod,a0)
        a1_now=ssp.copy_a(n_mod,a1)
        #
        # Fix maximum range!
        #
        statst=my.stats(pdl_flux)
        y_max_now=abs(statst[4]-statst[2])
        deft=1
        #
        # We forced a certain velocity!
        #
        if guided == 1:
            for j in range(0, n_mod):
                mask=gmask[ix,jy]
                #
                # FIXED velocity!
                #
                if link[j][3] == -1:
                    if typef[j] == "eline\n" and mask == 1:        
                        a_out_now[j][3]=gmap[ix,jy]
                        if guided_type == 0:
                            ia_now[j][3]=0
                        if guided_type == 1:
                            ia_now[j][3]=1
                            a0_now[j][3]=a_out_now[j][3]-30
                            a1_now[j][3]=a_out_now[j][3]+30
                #
                # FIXED DISPERSION
                #
                if link[j][2] == -1:
                    if guided_disp == 1:
                        if typef[j] == "eline\n":
                            val_disp=gmap_disp[ix,jy]
                            if val_disp > 0:
                                a_out_now[j][2]=val_disp
                                if guided_disp_type == 0:
                                    ia_now[j][2]=0
                                if guided_disp_type == 1:
                                    ia_now[j][2]=1
                                    a0_now[j][2]=a_out_now[j][2]-0.2*a_out_now[j][2]
                                    a1_now[j][2]=a_out_now[j][2]+0.2*a_out_now[j][2]
                if typef[j] == "eline\n":       
                    if deft == 1:
                        if link[j][1] == -1:
                            a1_max=2*y_max_now*(a[j][2]*((2*3.1416)**0.5));
                            a0_min=0.01*1.2*y_max_now*(a[j][2]*((2*3.1416)**0.5))
                            if a1[j][1] > a1_max:
                                a1_now[j][1]=a1_max
                                a0_now[j][1]=a0_min
        breakt=0
        if statst[0] != statst[1] or statst[0] != 0:
            n_max=5
            max_time=10
            i_max=0
            while i_ter < n_loops and breakt == 0:
                chi_iter=chi_sq_now
                pdl_grandom=ran.randn(nm)
                pdl_flux_fit=pdl_flux+pdl_grandom*pdl_e_flux*0.5
                [a_out_now, chi_sq_now,pdl_model]=ssp.fit_elines_grad_rnd_new_guided(pdl_wave,pdl_flux_fit,pdl_e_flux,n_mod,chi_goal,d_chi_goal,typef,a_out_now,ia_now,a0_now,a1_now,link,n_mc,pdl_masked,deft,scale_ini,guided_type,guided_disp_type)
                if chi_sq_now < chi_iter:
                ############################################################                                                       
                # Close to a result, narrow the range                                                                               
                    for i in range(0, n_mod):
                        for j in range(0, 9):
                            if typef[i] == "eline\n":
                                if ia_now[i][j] == 1:
                                    if link[i][j] == -1:
                                        a0_now_tmp=a0_now[i][j]
                                        a1_now_tmp=a1_now[i][j]
                                        delta_now=abs(a1_now[i][j]-a0_now[i][j])/(2.0)
                                        if j != 3:
                                            a0_now_tmp=a_out_now[i][j]-delta_now
                                            a1_now_tmp=a_out_now[i][j]+delta_now
                                        else:
                                            a0_now_tmp=a_out_now[i][j]-0.5*delta_now
                                            a1_now_tmp=a_out_now[i][j]+0.5*delta_now
                                        if a0_now_tmp < a0_now[i][j]:
                                            a0_now_tmp=a0_now[i][j]
                                        if a1_now_tmp > a1_now[i][j]:
                                            a1_now_tmp=a1_now[i][j]
                                        a0_now[i][j]=a0_now_tmp
                                        a1_now[i][j]=a1_now_tmp
                    ######################################                      
                    a_out_tmp=ssp.copy_a(n_mod,a_out_now)
                    a_results=ssp.copy_a_results(n_mod,nnk,a_out_tmp,a_results)    
                    nnk=nnk+1
                    i_ter=i_ter+1
                else:
                    #rnd_a=random[10]
                    i_max=i_max+1
                    if i_max > 5*n_loops:
                        breakt=1
            #plot=2
            if pdl_model.shape[0] == len(pdl_flux):
                pdl_model=np.transpose(pdl_model)
            pdl_model=pdl_model[0,:]
            if plot == 2:
                pdl_res=pdl_flux-pdl_model
                ssp.plot_results(2,pdl_wave,[pdl_flux,pdl_model,pdl_res],"junk","redshift")
        a_final=mean_a_results(n_mod,nnk,a_results,ix,jy,typef,pdl_wave)
        #print "- "+str(ix)+"/"+str(nx)+" "+str(jy)+"/"+str(ny)+" ---------------------------------------"
        #
        # Add background noise
        #
        pdl_res_now=pdl_flux-pdl_model
        stats_back=my.stats(pdl_res_now)
        a_final=ssp.add_back_noise(n_mod,a_final,typef,chi_sq_now,stats_back[1])
        #ssp.print_a_final(n_mod,a_final,typef,chi_sq_now)
        if jy == 0 and ix == 0:
            ssp.print_a_final_file(n_mod,a_final,typef,chi_sq_now,out_fit_spectra)
        else:
            ssp.print_a_final_file_add(n_mod,a_final,typef,chi_sq_now,out_fit_spectra,sile=1)
        org[:,ix,jy] = pdl_flux
        mod[:,ix,jy] = pdl_model
        res[:,ix,jy] = pdl_flux-pdl_model
        if memo == 1 and chi_sq_now < 2.0*chi_goal:
            a_out_now=ssp.copy_a_final(n_mod,a_final)    
        else:
            a_out_now=ssp.copy_a(n_mod,a_IN)
            ia_now=ssp.copy_a(n_mod,ia_IN)
            a0_now=ssp.copy_a(n_mod,a0_IN)
            a1_now=ssp.copy_a(n_mod,a1_IN)
h["CRPIX3"]=1
h["CRVAL3"]=crval+i0*cdelt
h["CDELT3"]=cdelt
h["NAXIS3"]=nm
wfits("kin_back_cube_org.fits", org, h)
wfits("kin_back_cube_res.fits", res, h)
wfits("kin_back_cube_mod.fits", mod, h)
#
# We create the maps!
#
print "Creating maps...";
for i in range(0, n_mod):
    if typef[i] != "poly1d\n":
        sec_f=pdl_f[i,:,:]
        sec_ef=pdl_ef[i,:,:]
        sec_v=pdl_v[i,:,:]
        sec_ev=pdl_ev[i,:,:]
        sec_dw=pdl_dw[i,:,:]
        sec_edw=pdl_edw[i,:,:]
        sec_dw_s=pdl_dw_s[i,:,:]
        sec_edw_s=pdl_edw_s[i,:,:]
        if i < 10:
            num="0"+str(i)
        else:
            num=str(i)
        outfile=prefix+"_"+num+".fits"
        outfile_f=prefix+"_flux_"+num+".fits"
        outfile_ef=prefix+"_eflux_"+num+".fits"
        outfile_w=prefix+"_vel_"+num+".fits"
        outfile_ew=prefix+"_evel_"+num+".fits"
        outfile_dw=prefix+"_disp_"+num+".fits"
        outfile_edw=prefix+"_edisp_"+num+".fits"
        outfile_dw_s=prefix+"_sigma_"+num+".fits"
        outfile_edw_s=prefix+"_esigma_"+num+".fits"
        wfits(outfile_f,sec_f,h)
        wfits(outfile_ef,sec_ef,h)
        wfits(outfile_w,sec_v,h)
        wfits(outfile_ew,sec_ev,h)
        wfits(outfile_dw,sec_dw,h)
        wfits(outfile_edw,sec_edw,h)
        wfits(outfile_dw_s,sec_dw_s,h)
        wfits(outfile_edw_s,sec_edw_s,h)
    else:
        sec_f=pdl_f[i,:,:]
        sec_ef=pdl_ef[i,:,:]
        outfile_f=prefix+"_cont.fits"
        outfile_ef=prefix+"_econt.fits"
        wfits(outfile_f,sec_f,h)
        wfits(outfile_ef,sec_ef,h)
sec_end=ssp.print_time()
sec_total=sec_end-sec_ini
print "# SECONDS = "+str(sec_total)
