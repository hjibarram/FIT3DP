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
import math

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
if len(sys.argv) < 3:
    print "USE: map_auto_ssp_rnd_seg.py elines_OUT seg_file.fits PREFIX_OUT [DISP_INST] [WAVE_NORM]"
    sys.exit(0)

infile=sys.argv[1]
sspfile=sys.argv[1]
cut="elines_"
sspfile=sspfile.replace(cut,'')
seg_file=sys.argv[2]
prefix=sys.argv[3]

disp0=2.45
wave0=4275
if len(sys.argv) == 5:
    disp0=float_(sys.argv[4])
if len(sys.argv) == 6:
    disp0=float_(sys.argv[4])
    wave0=float_(sys.argv[5])

[pdl_seg, hdr]=gdata(seg_file, 0, header=True)
[nx,ny]=pdl_seg.shape

call="map_auto_ssp_AGE_MET_rnd_seg.py "+sspfile+" "+seg_file+" "+prefix+" "+str(wave0)+" "+str(disp0)
mycall(call)

print "Reading fit output file"
n=0
l=0
nmod=0
N=0
first=0.0;
flux=np.zeros([1,1])
eflux=np.zeros([1,1])
sig=np.zeros([1,1])
esig=np.zeros([1,1])
vel=np.zeros([1,1])
evel=np.zeros([1,1])
wave=[]
wave_name=[]
reg_na=[]
stp=0
f=open(infile,"r")
for line in f:
    data=line.split(" ")
    data=filter(None,data)
    if not "#" in data[0]:
        if "eline" in data[0]:
            if first > 0 and first == float_(data[1]):
                n=n+1
                nmod=l
                flux=np.concatenate((flux,np.zeros([1,l])),axis=0)
                eflux=np.concatenate((eflux,np.zeros([1,l])),axis=0)
                sig=np.concatenate((sig,np.zeros([1,l])),axis=0)
                esig=np.concatenate((esig,np.zeros([1,l])),axis=0)
                vel=np.concatenate((vel,np.zeros([1,l])),axis=0)
                evel=np.concatenate((evel,np.zeros([1,l])),axis=0)
                l=0
                stp=1
            if n == 0 and l==0:
                first=float_(data[1])
            if l > 0 and n == 0 and stp == 0:
                flux=np.concatenate((flux,np.zeros([1,1])),axis=1)
                eflux=np.concatenate((eflux,np.zeros([1,1])),axis=1)
                sig=np.concatenate((sig,np.zeros([1,1])),axis=1)
                esig=np.concatenate((esig,np.zeros([1,1])),axis=1)
                vel=np.concatenate((vel,np.zeros([1,1])),axis=1)
                evel=np.concatenate((evel,np.zeros([1,1])),axis=1)
            if n == 0:
                wave.extend([float_(data[1])])
                [dw,w]=math.modf(wave[l])
                wave_name.extend([str(w)])
            else:
                wave[l]=float_(data[1])
                [dw,w]=math.modf(wave[l])
                wave_name[l]=str(w)
            flux[n][l]=float_(data[3])
            eflux[n][l]=float_(data[4])
            sig[n][l]=float_(data[5])
            esig[n][l]=float_(data[6])
            vel[n][l]=float_(data[7])
            evel[n][l]=float_(data[8])
            l=l+1
    if "#ID" in line:
        reg_na.extend([float_(line.replace("\n","").split(" ")[1])])
reg_na=np.array(reg_na)
f.close()
print str(nmod)+", "+str(n)+" data"
print "DONE"
print "Feeding arrays"
pdl_flux=np.zeros([nmod,nx,ny])
pdl_eflux=np.zeros([nmod,nx,ny])
pdl_vel=np.zeros([nmod,nx,ny])
pdl_disp=np.zeros([nmod,nx,ny])
pdl_disp_km_s=np.zeros([nmod,nx,ny])
N=n+1
for ix in range(0, nx):
    for iy in range(0, ny):
        iseg=int(pdl_seg[ix,iy])
        if iseg > 0:
            iseg=iseg-1
            for i in range(0, nmod):
                if iseg in reg_na:
                    iseg2=np.where(reg_na == iseg)[0][0]
                    pdl_flux[i,ix,iy]=flux[iseg2][i]
                    pdl_eflux[i,ix,iy]=eflux[iseg2][i]
                    pdl_disp[i,ix,iy]=sig[iseg2][i]
                    if sig[iseg2][i] > disp0:
                        sig_km_s=(np.sqrt(sig[iseg2][i]**2.0-disp0**2.0)/float_(wave_name[i]))*300000
                    else:
                        sim_km_s=0
                    pdl_disp_km_s[i,ix,iy]=sig_km_s
                    pdl_vel[i,ix,iy]=vel[iseg2][i]
print "Done"
print "Writing output"
for i in range(0, nmod):
    ii=i
    if ii < 10:
        iis="0"+str(i)
    flux_name=prefix+"_flux_"+wave_name[i]+".fits"
    eflux_name=prefix+"_eflux_"+wave_name[i]+".fits"
    disp_name=prefix+"_disp_"+wave_name[i]+".fits"
    disp_km_s_name=prefix+"_disp_km_h_"+wave_name[i]+".fits"
    vel_name=prefix+"_vel_"+wave_name[i]+".fits"
    map_flux=pdl_flux[i,:,:]
    map_eflux=pdl_eflux[i,:,:]
    map_vel=pdl_vel[i,:,:]
    map_disp=pdl_disp[i,:,:]
    map_disp_km_s=pdl_disp_km_s[i,:,:]
    wfits(flux_name,map_flux,hdr)
    wfits(eflux_name,map_eflux,hdr)
    wfits(vel_name,map_vel,hdr)
    wfits(disp_name,map_disp,hdr)
    wfits(disp_km_s_name,map_disp_km_s,hdr)
    call="write_img_header.py "+flux_name+" CRVAL1 "+str(wave[i])
    mycall(call)
    call="write_img_header.py "+eflux_name+" CRVAL1 "+str(wave[i])
    mycall(call)
    call="write_img_header.py "+vel_name+" CRVAL1 "+str(wave[i])
    mycall(call)
    call="write_img_header.py "+disp_name+" CRVAL1 "+str(wave[i])
    mycall(call)
    call="write_img_header.py "+disp_km_s_name+" CRVAL1 "+str(wave[i])
    mycall(call)
    print str(wave[i])+" DONE"
print "DONE"