#! /usr/bin/env python
import sys
import numpy as np
from numpy import float_
from numpy import absolute as abs
from numpy import random as ran
import cosmolopy.distance as cd
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
    import os
    os.system(comand)

def wfits(name, data, hdr):
    if ptt.exists(name) == False:
        wfit(name,data,hdr)
    else:
        sycall("rm "+name)
        wfit(name,data,hdr)
        
sys.argv=filter(None,sys.argv)
if len(sys.argv) < 2:
    print "USE: sum_mass_age.py NAME"
    sys.exit(0)
name=sys.argv[1]
seg="cont_seg."+name+".fits.gz"
pdl_seg=gdata(seg)
[nx,ny]=pdl_seg.shape
flux_file="map.CS."+name+"_flux_ssp.fits.gz"
[pdl_flux, hdr]=gdata(flux_file, 0, header=True)
[nx,ny]=pdl_flux.shape
maskName="mask."+name+".V.fits.gz"
#maskName="DMASK."+name+".fits.gz"
if ptt.exists(maskName) == True:
    m_file=maskName
    #pdl_mask=1-gdata(m_file)
    pdl_mask=gdata(m_file)
else:
    pdl_mask=np.ones([nx,ny])
nl=0
med_flux=[]
e_med_flux=[]
a_redshift=[]
sum_flux=0
file="auto_ssp.CS."+name+".rss.out"
f=open(file, "r")
for line in f:
    if not "#" in line:
        data=line.split(",")
        data=filter(None,data)
        med_flux.extend([float_(data[13])])
        e_med_flux.extend([float_(data[14])])
        a_redshift.extend([float_(data[7])])
        sum_flux=sum_flux+float_(data[13])
        nl=nl+1
f.close()
med_flux=np.array(med_flux)
e_med_flux=np.array(e_med_flux)
a_redshift=np.array(a_redshift)
sum_flux=sum_flux/float_(nl)
redshift=np.median(a_redshift)
cosmo = {'omega_M_0' : 0.27, 'omega_lambda_0' : 0.73, 'h' : 0.71}
cosmo = cd.set_omega_k_0(cosmo)
DL=cd.luminosity_distance(redshift,**cosmo)
ratio=3.08567758e24
modz=5.0*np.log10(DL)+25.0
DL=DL*ratio
L=4.0*np.pi*(DL**2.0)#/(1+$redshift);
Factor=(L*1e-16)/3.826e33
file="coeffs_auto_ssp.CS."+name+".rss.out"
out_pack_csv="sum_mass_age."+name+".pack.csv"
f=open(out_pack_csv,"w")
f.write("# (1) ID\n")
f.write("# (2) File\n")
f.write("# (3) Description\n")
f.write("# (4) Type\n")
f.write("# (5) units\n")
a_age=[]
a_met=[]
f2=open(file,"r")
n=0
n_ini=0
n_age=0
n_met=0
ID=np.zeros(156)
AGE=np.zeros(156)
MET=np.zeros(156)
COEFF=np.zeros(156)
NORM=np.zeros(156)
ML=np.zeros(156)
AV=np.zeros(156)
SUM=np.zeros(156)
SUM_MASS=np.zeros(156)
for line in f2:
    if n_ini < 156:
        if not "#" in line:
            data=line.split(" ")
            data=filter(None,data)
            n=int(data[0])
            ID[n]=int(data[0])
            AGE[n]=float_(data[1])
            MET[n]=float_(data[2])
            COEFF[n]=float_(data[3])
            NORM[n]=float_(data[4])
            ML[n]=float_(data[5])
            AV[n]=float_(data[6])
            SUM[n]=SUM[n]+float_(data[3])
            SUM_MASS[n]=SUM_MASS[n]+10.0**((-1)*float_(data[5]))
            diff_age=1
            for i in range(0, n):
                if AGE[n] == AGE[i]:
                    diff_age=0
            if diff_age == 1:
                a_age.extend([AGE[n]])
                n_age=n_age+1
            diff_met=1
            for i in range(0, n):
                if MET[n] == MET[i]:
                    diff_met=0
            if diff_met == 1:
                a_met.extend([MET[n]])
                n_met=n_met+1
            n_ini=n_ini+1
f2.close()
n=n+1
print "# N="+str(n_ini)
print "# N.AGE="+str(n_age)
print "# N.MET="+str(n_met)
pdl_ML=np.zeros([nx,ny])
pdl_Mass=np.zeros([nx,ny])
pdl_e_Mass=np.zeros([nx,ny])
i_n=1


for i in range(0, n):
    norm_file="map.CS."+name+"_NORM_"+str(i)+"_ssp.fits.gz"
    pdl_tmp=gdata(norm_file)
    norm_e_file="map.CS."+name+"_eNORM_"+str(i)+"_ssp.fits.gz"
    pdl_e_tmp=gdata(norm_e_file)
    pdl_ML=pdl_ML+pdl_tmp*10.0**(ML[i])
    pdl_Mass=pdl_Mass+pdl_tmp*10.0**(ML[i])*pdl_flux# /$mstar[$i];
    pdl_e_Mass=pdl_e_Mass+pdl_e_tmp*10.0**(ML[i])*pdl_flux#/$mstar[$i];
    if i == 0:
        [nx,ny]=pdl_tmp.shape
        pdl_cube=np.zeros([n,nx,ny])
        pdl_e_cube=np.zeros([n,nx,ny])
    pdl_cube[i,:,:] = pdl_tmp
    pdl_e_cube[i,:,:] = pdl_e_tmp
    file_name="map.CS.NAME_NORM_"+str(i)+"_ssp.fits.gz"
    f.write(str(i_n)+","+file_name+",Luminosity Fraction for age-met "+str(AGE[i])+"-"+str(MET[i])+" SSP, flux, fraction\n")
    i_n=i_n+1
pdl_cube=pdl_cube*pdl_mask
pdl_e_cube=pdl_e_cube*pdl_mask
file_Av="map.CS."+name+"_Av_ssp.fits.gz"
pdl_Av=gdata(file_Av)
pdl_Mass=pdl_Mass*Factor#+$pdl_tmp*10**($ML[$i])*$pdl_flux/$mstar[$i];
pdl_e_Mass=pdl_e_Mass*Factor#+$pdl_tmp*10**($ML[$i])*$pdl_flux/$mstar[$i];
pdl_Mass[np.isnan(pdl_Mass)]=1
pdl_Mass[np.where(pdl_Mass == 0)]=1
pdl_l_Mass=np.log10(pdl_Mass)
pdl_l_Mass=pdl_l_Mass*pdl_mask
mass_file="map.CS."+name+"_Mass_ssp.fits"
wfits(mass_file, pdl_l_Mass , hdr)
call="gzip -f "+mass_file
sycall(call)
pdl_e_l_Mass=0.4*pdl_e_Mass/pdl_Mass;
pdl_e_l_Mass=pdl_e_l_Mass*pdl_mask
mass_e_file="map.CS."+name+"_eMass_ssp.fits"
wfits(mass_e_file, pdl_e_l_Mass , hdr)
call="gzip -f "+mass_e_file
sycall(call)
mass_file="map.CS."+name+"_Mass_dust_cor_ssp.fits"
pdl_l_Mass_dust_corr=pdl_l_Mass+0.4*pdl_Av
wfits(mass_file, pdl_l_Mass_dust_corr, hdr)
call="gzip -f "+mass_file
sycall(call)
mass_file="map.CS."+name+"_ML_ssp.fits"
pdl_ML[np.isnan(pdl_ML)]=1
pdl_ML[np.where(pdl_ML == 0)]=1
pdl_l_ML=np.log10(pdl_ML)
wfits(mass_file, pdl_l_ML, hdr)
call="gzip -f "+mass_file
sycall(call)
Mass_sum=np.sum(pdl_Mass)
e_Mass_sum=np.sum(pdl_e_Mass)
l_Mass=np.log10(Mass_sum)
e_l_Mass=0.4*e_Mass_sum/Mass_sum
Mass_sum_dust_corr=0
pdl_l_Mass_dust_corr[np.isnan(pdl_l_Mass_dust_corr)]=0
pdl_l_Mass[np.isnan(pdl_l_Mass)]=0
for i in range(0, nx):
    for j in range(0, ny):
        val=abs(pdl_l_Mass_dust_corr[i,j])
        val2=abs(pdl_l_Mass[i,j])
        if val > 0:
            Mass_sum_dust_corr=Mass_sum_dust_corr+10.0**(val)
        else:
            if val2 > 0:
                Mass_sum_dust_corr=Mass_sum_dust_corr+10.0**(val2)
l_Mass_dust_corr=np.log10(Mass_sum_dust_corr)
print "# name, log10(Mass)"
print "Mass,"+name+","+str(l_Mass)+","+str(l_Mass_dust_corr)+","+str(e_l_Mass)
for i in range(0, n_age):
    age_now=a_age[i]
    age_t=age_now
    if age_t < 10:
        age_t="0"+("%6.4f" % age_t)
    else:
        age_t=("%6.4f" % age_t)
    age_file="map.CS."+name+"_"+age_t+"_NORM_age.fits"
    age_e_file="map.CS."+name+"_"+age_t+"_eNORM_age.fits"
    pdl_age=np.zeros([nx,ny])
    pdl_e_age=np.zeros([nx,ny])
    for j in range(0, n):
        if age_now == AGE[j]:
            pdl_age=pdl_age+pdl_cube[j,:,:]
            pdl_e_age=pdl_e_age+pdl_e_cube[j,:,:]**2.0
    wfits(age_file, pdl_age, hdr)
    pdl_e_age=np.sqrt(pdl_e_age)
    wfits(age_e_file, pdl_e_age, hdr)
    call="gzip -f "+age_file
    sycall(call)
    call="gzip -f "+age_e_file
    sycall(call)
    file_name="map.CS.NAME_"+age_t+"_NORM_age.fits.gz"
    f.write(str(i_n)+","+file_name+",Luminosity Fraction for age "+age_t+" SSP, flux, fraction\n")
    i_n=i_n+1
for i in range(0, n_met):
    met_now=a_met[i]
    met_file="map.CS."+name+"_"+("%6.4f" % met_now)+"_NORM_met.fits"
    pdl_met=np.zeros([nx,ny])
    met_e_file="map.CS."+name+"_"+("%6.4f" % met_now)+"_eNORM_met.fits"
    pdl_e_met=np.zeros([nx,ny])
    for j in range(0, n):
        if met_now == MET[j]:
            pdl_met=pdl_met+pdl_cube[j,:,:]
            pdl_e_met=pdl_e_met+pdl_e_cube[j,:,:]**2.0
    wfits(met_file, pdl_met, hdr)
    call="gzip -f "+met_file
    sycall(call)
    pdl_e_met=np.sqrt(pdl_e_met)
    wfits(met_e_file, pdl_e_met, hdr)
    call="gzip -f "+met_e_file
    sycall(call)
    file_name="map.CS.NAME_"+("%6.4f" % met_now)+"_NORM_met.fits.gz"
    f.write(str(i_n)+","+file_name+",Luminosity Fraction for met "+("%6.4f" % MET[i])+" SSP, flux, fraction\n")
    i_n=i_n+1
f.close()