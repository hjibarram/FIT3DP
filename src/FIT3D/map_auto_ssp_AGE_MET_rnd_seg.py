#!/usr/bin/python
import sys
import numpy as np
from numpy import float_
from numpy import absolute as abs
import pyfits as pyf
from pyfits import getheader as ghead
from pyfits import getdata as gdata
from pyfits import writeto as wfit
import os.path as ptt
import my
from my import mycall

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

vel_light=299792.458
sys.argv=filter(None,sys.argv)
if len(sys.argv) < 6:
    print "USE: map_auto_ssp.py AUTO.OUT seg_file.fits  PREFIX_OUT WAVE_CEN INST_DIST"
    sys.exit(0)
infile=sys.argv[1]
seg_file=sys.argv[2]
prefix=sys.argv[3]
wave_cen=float_(sys.argv[4])
inst_disp=float_(sys.argv[5])
[pdl_seg, hdr]=gdata(seg_file, 0, header=True)
[nx,ny]=pdl_seg.shape
pdl_chi=np.zeros([nx,ny])
pdl_age=np.zeros([nx,ny])
pdl_age_mass=np.zeros([nx,ny])
pdl_e_age=np.zeros([nx,ny])
pdl_age_lum=np.zeros([nx,ny])
pdl_e_age_lum=np.zeros([nx,ny])
pdl_age_lum_M=np.zeros([nx,ny])
pdl_e_age_lum_M=np.zeros([nx,ny])
pdl_met_lum=np.zeros([nx,ny])
pdl_met_mass=np.zeros([nx,ny])
pdl_met_lumI=np.zeros([nx,ny])
pdl_met_lumII=np.zeros([nx,ny])
pdl_met_lumII_mass=np.zeros([nx,ny])
pdl_e_met_lumII=np.zeros([nx,ny])
pdl_met_e_lumII_mass=np.zeros([nx,ny])
pdl_met=np.zeros([nx,ny])
pdl_e_met=np.zeros([nx,ny])
pdl_Av=np.zeros([nx,ny])
pdl_e_Av=np.zeros([nx,ny])
pdl_disp=np.zeros([nx,ny])
pdl_vel=np.zeros([nx,ny])
pdl_e_disp=np.zeros([nx,ny])
pdl_e_vel=np.zeros([nx,ny])
pdl_flux=np.zeros([nx,ny])
pdl_e_flux=np.zeros([nx,ny])
pdl_ml=np.zeros([nx,ny])
pdl_mass=np.zeros([nx,ny])
pdl_disp_km=np.zeros([nx,ny])
pdl_e_disp_km=np.zeros([nx,ny])
C_lim=1/100.
print "Reading input files"
nc=-1
nc_max=-1
SUM=[]
SUM_MASS=[]
SUM_now=[]
AGE=np.zeros([1,1])
MET=np.zeros([1,1])
COEFF=np.zeros([1,1])
NORM=np.zeros([1,1])
eNORM=np.zeros([1,1])
ML=np.zeros([1,1])
AV=np.zeros([1,1])
COEFF_MASS=np.zeros([1,1])
coeffs_infile="coeffs_"+infile
f=open(coeffs_infile,"r")
for line in f:
    data=line.split(" ")
    data=filter(None,data)
    if data[0] != "#":
        i=int(data[0])
        if i == 0:
            nc=nc+1
            SUM.extend([0])
            SUM_MASS.extend([0])
            SUM_now.extend([0])
            if nc > 0:
                if nc == 1:
                    lsp=AGE.shape[0]
                AGE=np.concatenate((AGE,np.zeros([lsp,1])),axis=1)
                MET=np.concatenate((MET,np.zeros([lsp,1])),axis=1)
                COEFF=np.concatenate((COEFF,np.zeros([lsp,1])),axis=1)
                NORM=np.concatenate((NORM,np.zeros([lsp,1])),axis=1)
                eNORM=np.concatenate((eNORM,np.zeros([lsp,1])),axis=1)
                ML=np.concatenate((ML,np.zeros([lsp,1])),axis=1)
                AV=np.concatenate((AV,np.zeros([lsp,1])),axis=1)
                COEFF_MASS=np.concatenate((COEFF_MASS,np.zeros([lsp,1])),axis=1)                
        if nc == 0 and i > 0:
            AGE=np.concatenate((AGE,np.zeros([1,1])),axis=0)
            MET=np.concatenate((MET,np.zeros([1,1])),axis=0)
            COEFF=np.concatenate((COEFF,np.zeros([1,1])),axis=0)
            NORM=np.concatenate((NORM,np.zeros([1,1])),axis=0)
            eNORM=np.concatenate((eNORM,np.zeros([1,1])),axis=0)
            ML=np.concatenate((ML,np.zeros([1,1])),axis=0)
            AV=np.concatenate((AV,np.zeros([1,1])),axis=0)
            COEFF_MASS=np.concatenate((COEFF_MASS,np.zeros([1,1])),axis=0)
        AGE[i][nc]=float_(data[1])
        MET[i][nc]=float_(data[2])
        COEFF[i][nc]=float_(data[3])
        NORM[i][nc]=float_(data[4])
        eNORM[i][nc]=float_(data[8])
        ML[i][nc]=float_(data[5])
        AV[i][nc]=float_(data[6])
        SUM[nc]=SUM[nc]+float_(data[3])
        COEFF_MASS[i][nc]=(float_(data[3])*10**(float_(data[5])))
        SUM_MASS[nc]=SUM_MASS[nc]+(float_(data[3])*10**(float_(data[5])))
        if COEFF[i][nc] > C_lim:
            SUM_now[nc]=SUM_now[nc]+COEFF[i][nc]
        if nc_max <= i:
            nc_max=i+1
f.close()
pdl_COEFF=np.zeros([nc_max,nx,ny])
pdl_NORM=np.zeros([nc_max,nx,ny])    
pdl_eNORM=np.zeros([nc_max,nx,ny])  
pdl_AV=np.zeros([nc_max,nx,ny])   
print "DONE"
chi=[]
age=[]
e_age=[]
met=[]
e_met=[]
Av=[]
e_Av=[]
vel=[]
e_vel=[]
disp=[]
e_disp=[]
flux=[]
e_flux=[]
n=0
f=open(infile, "r")
for line in f:
    data=line.split(",")
    data=filter(None,data)
    if not "#" in data[0]:
        chi.extend([float_(data[0])])
        age.extend([float_(data[1])])
        e_age.extend([float_(data[2])])
        met.extend([float_(data[3])])
        e_met.extend([float_(data[4])])
        Av.extend([float_(data[5])])
        e_Av.extend([float_(data[6])])
        vel.extend([float_(data[7])*299792.458])
        e_vel.extend([float_(data[8])*299792.458])
        disp.extend([float_(data[9])/(1+float_(data[7]))])
        e_disp.extend([abs(float_(data[10]))])
        flux.extend([float_(data[13])])
        e_flux.extend([float_(data[14])])
        n=n+1
f.close()
print "Feeding the arrays"
for ix in range(0, nx):
    for iy in range(0, ny):
        iseg=int(pdl_seg[ix,iy])
        if iseg > 0:
            iseg=iseg-1
            pdl_chi[ix,iy]=chi[iseg]
            pdl_age[ix,iy]=age[iseg]
            pdl_met[ix,iy]=met[iseg]
            pdl_Av[ix,iy]=Av[iseg]
            pdl_vel[ix,iy]=vel[iseg]
            pdl_disp[ix,iy]=disp[iseg]
            if disp[iseg] > inst_disp:
                disp_km=(np.sqrt(abs(disp[iseg]**2.-inst_disp**2.))/wave_cen)*299792.458
            else:
                disp_km=0
            pdl_disp_km[ix,iy]=disp_km
            pdl_flux[ix,iy]=flux[iseg]
            pdl_e_age[ix,iy]=e_age[iseg]
            pdl_e_met[ix,iy]=e_met[iseg]
            pdl_e_Av[ix,iy]=e_Av[iseg]
            pdl_e_vel[ix,iy]=e_vel[iseg]
            pdl_e_disp[ix,iy]=e_disp[iseg]
            e_disp_km=(e_disp[iseg]/wave_cen)*299792.458
            pdl_e_disp_km[ix,iy]=e_disp_km
            pdl_e_flux[ix,iy]=e_flux[iseg]
            for i in range(0, nc_max):
                pdl_NORM[i,ix,iy]=COEFF[i][iseg]
                pdl_eNORM[i,ix,iy]=eNORM[i][iseg]
                age_now=AGE[i][0]
                met_now=MET[i][0]
                age_val=pdl_age_lum[ix,iy]        
                age_val=age_val+COEFF[i][iseg]*(np.log10(age_now))
                pdl_age_lum[ix,iy]=age_val
                if SUM_MASS[iseg] > 0:
                    age_val=pdl_age_mass[ix,iy]        
                    age_val=age_val+COEFF_MASS[i][iseg]*(np.log10(age_now))/SUM_MASS[iseg]
                    pdl_age_mass[ix,iy]=age_val
                    ml_now=pdl_ml[ix,iy]
                    ml_now=ml_now+COEFF[i][iseg]*ML[i][iseg]
                    pdl_ml[ix,iy]=ml_now
                    mass_now=pdl_mass[ix,iy]
                    mass_now=mass_now+COEFF[i][iseg]*10**(ML[i][iseg])*flux[iseg]*1e-16
                    pdl_mass[ix,iy]=mass_now
                met_val=pdl_met_lum[ix,iy]        
                met_val=met_val+COEFF[i][iseg]*met_now
                pdl_met_lum[ix,iy]=met_val
                met_val=pdl_met_lumI[ix,iy]     
                met_val=met_val+NORM[i][iseg]*np.log10(met_now/0.02)
                pdl_met_lumI[ix,iy]=met_val
                met_val=pdl_met_lumII[ix,iy]       
                met_val=met_val+COEFF[i][iseg]*np.log10(met_now/0.02)
                pdl_met_lumII[ix,iy]=met_val
                if SUM_MASS[iseg] > 0:
                    met_val=pdl_met_mass[ix,iy]                
                    met_val=met_val+COEFF_MASS[i][iseg]*np.log10(met_now/0.02)/SUM_MASS[iseg]
                    pdl_met_mass[ix,iy]=met_val
            met_now=pdl_met[ix,iy]
            met_e_now=pdl_e_met[ix,iy]
            if met_now > 0:
                val_m=met_e_now/met_now
                pdl_e_met_lumII[ix,iy]=val_m
            age_now=pdl_age_lum[ix,iy]
            e_age_now=pdl_e_age[ix,iy]
            if age_now > 0:
                e_age_now=e_age_now/age_now
            pdl_e_age_lum[ix,iy]=age_val
            age_val=9+pdl_age_lum[ix,iy]
            pdl_age_lum[ix,iy]=age_val
            age_val=9+pdl_age_mass[ix,iy]
            pdl_age_lum_M[ix,iy]=age_val
            met_val=pdl_met_lum[ix,iy]        
            met_val=np.log10(met_val/0.02)
            pdl_met_lum[ix,iy]=met_val
print "DONE"
print "Writting the output files"
chi_name=prefix+"_chi_ssp.fits"
wfits(chi_name,pdl_chi,hdr)
age_name=prefix+"_age_ssp.fits"
wfits(age_name,pdl_age,hdr)
age_name=prefix+"_age_mass_ssp.fits"
wfits(age_name,pdl_age_mass,hdr)
e_age_name=prefix+"_e_age_ssp.fits"
wfits(e_age_name,pdl_e_age,hdr)
age_name=prefix+"_log_age_yr_ssp.fits"
wfits(age_name,pdl_age_lum,hdr)
age_name=prefix+"_e_log_age_yr_ssp.fits"
wfits(age_name,pdl_e_age_lum,hdr)
age_name=prefix+"_log_age_yr_mass_ssp.fits"
wfits(age_name,pdl_age_lum_M,hdr)
age_name=prefix+"_e_log_age_yr_ssp.fits"
wfits(age_name,pdl_e_age_lum,hdr)
met_name=prefix+"_met_ssp.fits"
wfits(met_name,pdl_met,hdr)
met_name=prefix+"_met_ZH_mass_ssp.fits"
wfits(met_name,pdl_met_mass,hdr)
met_name=prefix+"_met_ZH_ssp.fits"
wfits(met_name,pdl_met_lumII,hdr)
met_name=prefix+"_e_met_ZH_ssp.fits"
wfits(met_name,pdl_e_met_lumII,hdr)
Av_name=prefix+"_Av_ssp.fits"
wfits(Av_name,pdl_Av,hdr)
disp_name=prefix+"_disp_ssp.fits"
wfits(disp_name,pdl_disp,hdr)
disp_km_name=prefix+"_disp_km_h_ssp.fits"
wfits(disp_km_name,pdl_disp_km,hdr)
vel_name=prefix+"_vel_ssp.fits"
wfits(vel_name,pdl_vel,hdr)
flux_name=prefix+"_flux_ssp.fits"
wfits(flux_name,pdl_flux,hdr)
e_met_name=prefix+"_e_met_ssp.fits"
wfits(e_met_name,pdl_e_met,hdr)
e_Av_name=prefix+"_e_Av_ssp.fits"
wfits(e_Av_name,pdl_e_Av,hdr)
e_disp_name=prefix+"_e_disp_ssp.fits"
wfits(e_disp_name,pdl_e_disp,hdr)
e_disp_km_name=prefix+"_e_disp_km_h_ssp.fits"
wfits(e_disp_km_name,pdl_e_disp_km,hdr)
e_vel_name=prefix+"_e_vel_ssp.fits"
wfits(e_vel_name,pdl_e_vel,hdr)
e_flux_name=prefix+"_e_flux_ssp.fits"
wfits(e_flux_name,pdl_e_flux,hdr)
mass_name=prefix+"_mass_ssp.fits"
wfits(mass_name,pdl_mass,hdr)
ml_name=prefix+"_ml_ssp.fits"
wfits(ml_name,pdl_ml,hdr)

for j in range(0, nc_max):
    file_name=prefix+"_NORM_"+str(j)+"_ssp.fits"
    sec=pdl_NORM[j,:,:]
    wfits(file_name,sec,hdr)
    call="write_img_header.py "+file_name+" AGE "+str(AGE[j][0])
    mycall(call)
    call="write_img_header.py "+file_name+" MET "+str(MET[j][0])
    mycall(call)
    
    file_name=prefix+"_eNORM_"+str(j)+"_ssp.fits"
    sec=pdl_eNORM[j,:,:]
    wfits(file_name,sec,hdr)
    call="write_img_header.py "+file_name+" AGE "+str(AGE[j][0])
    mycall(call)
    call="write_img_header.py "+file_name+" MET "+str(MET[j][0])
    mycall(call)