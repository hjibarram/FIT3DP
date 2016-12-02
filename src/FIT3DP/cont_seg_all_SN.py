#! /usr/bin/env python
import sys
import numpy as np
from numpy import float_
import pyfits as py
from pyfits import getheader as ghead
from pyfits import getdata as gdata
from pyfits import writeto as wfits
import os.path as pt
import my

def sycall(comand):
    from subprocess import call
    line=comand.split(" ")
    fcomand=[]
    fcomand.extend(line)
    call(fcomand)

def stats(data):
    out=np.zeros(7)
    out[0]=np.mean(data)
    out[1]=np.mean(data)+np.std(data)
    out[2]=np.median(data)
    out[3]=np.amin(data)
    out[4]=np.amax(data)
    out[5]=np.std(data)
    out[6]=np.mean(data)+np.std(data)
    return out

def is_max(i,j):
    i1=i-1
    i2=i+2
    j1=j-1
    j2=j+2
    if i1 < 0:
        i1=0
        i2=i1+3
    if i2 > nx-1:
        i2=nx-1
        i1=i2-3
    if j1 < 0:
        j1=0
        j2=j1+3
    if j2 > ny-1:
        j2=ny-1
        j1=j2-3
    a=[]
    k=0
    is_max=-3
    for ii in range(i1, i2):
        for jj in range(j1, j2):
            if ii !=i and jj != j:
                if ii > 0 and ii < nx and jj > 0 and jj < ny:
                    nr=seg_map[ii,jj]
                    if nr == 0:
                        a.extend([Ha_map[ii,jj]])
                        k=k+1
    if len(a) == 0:
        a=[0,0]
    a=np.array(a)
    med_a=np.median(a)
    sig_a=np.std(a)
    min=np.amin(a)
    max=np.amax(a)
    if i > 0 and i < nx and j > 0 and j < ny:
        val=Ha_map[i,j]
        F_peak=Ha_map[ip,jp]
        if val > frac*F_peak:
            ncut=2.0
            if val-med_a < 0:
                is_max=-1
            if abs(med_a-val) <= ncut*sig_a:
                is_max=0
            if val > med_a:
                is_max=1
    return is_max

def is_near(i,j,seg_map,nx,ny):
    i1=i-1
    i2=i+2
    j1=j-1
    j2=j+2
    if i1 < 0:
        i1=0
        i2=i1+3
    if i2 > nx-1:
        i2=nx-1
        i1=i2-3
    if j1 < 0:
        j1=0
        j2=j1+3
    if j2 > ny-1:
        j2=ny-1
        j1=j2-3
    a=0
    k=0
    is_max=-3
    for ii in range(i1, i2):
        for jj in range(j1, j2):
            if ii != i and jj != j:
                if ii > 0 and ii < nx and jj > 0 and jj < ny:
                    val=seg_map[ii,jj]
                    if val > 0:
                        a=a+1
                k=k+1
    R=a
    return R

def DO(i,j,mask_map,seg_map):
    dist=np.sqrt((i-ip)**2+(j-jp)**2)
    #is_maxt=is_max(i,j,Ha_map,seg_map)
    F_PEAK_NOW=Ha_map[ip,jp]
    if NR == 1:
        frac_dist=1
    else:
        if (F_PEAK_NOW/F_PEAK[1])**0.5 !=0:
            frac_dist=0.5/((F_PEAK_NOW/F_PEAK[1])**0.5)
        else:
            frac_dist=np.nan
        if frac_dist > 1:
            frac_dist=1
    dist_max_now=max_dist/frac_dist
    if dist_max_now < 1.5:
        dist_max_now=1.5
    is_point=0
    if i > 0 and i < nx and j > 0 and j < ny:
        val=Ha_map[i,j]
        F_peak=Ha_map[ip,jp]
        if dist < max_dist and val > MIN_FLUX:
            if dist < max_dist and val > frac*F_peak:
                nr=seg_map[i,j]
                if nr == 0:
                    seg_map[i,j]=NR
                    is_point=1 
                mask_map[i,j]=0        
    return [is_point,mask_map,seg_map]

def LOOP(ist,js,seg_map,mask_map):
    nmax=(nx-1)*(ny-1)
    n_now=0
    delta=2
    s=0
    k=0
    sig=0
    i=ist
    j=js
    seg_map[i,j]=NR
    mask_map[i,j]=0
    n_now=n_now+s
    breakt=0
    while breakt == 0:
        if sig == 0:
            for j in range(js+1, js+delta):
                [s,mask_map,seg_map]=DO(i,j,mask_map,seg_map)
                n_now=n_now+1
            j=j-1
            for i in range(ist+1, ist+delta):
                [s,mask_map,seg_map]=DO(i,j,mask_map,seg_map)
                n_now=n_now+1
            sig=1
            i=i-1
        else:
            for j in range(js-1, js-delta,-1):
                [s,mask_map,seg_map]=DO(i,j,mask_map,seg_map)
                n_now=n_now+1
            j=j+1
            for i in range(ist-1, ist-delta,-1):
                [s,mask_map,seg_map]=DO(i,j,mask_map,seg_map)
                n_now=n_now+1
            sig=0
            i=i+1
        delta=delta+1
        js=j
        ist=i
        dist=np.sqrt((i-ip)**2+(j-jp)**2)
        if i > 0 and j > 0 and i < nx and j < ny:
            val_now=map_now[i,j]
        else:
            val_now=0
        if dist > 1.5*max_dist or val_now > 10*F_max:
            breakt=1
    return [seg_map,mask_map]

def nearest(i,j,NR,IP,JP,seg_map):
    min_dist=1e12
    for ii in range(0, NR):
        dist=np.sqrt((IP[ii]-i)**2+(JP[ii]-j)**2)
        if dist < min_dist:
            min_dist=dist
            seg_map[i,j]=ii
    return seg_map

sys.argv=filter(None,sys.argv)
if (sys.argv) < 8:
    print "USE: cont_seg_all_SN.py signal.fits,noise.fits FLUX_LIMIT_PEAK TARGET_SN,MIN_SN FRAC_PEAK MIN_FLUX SEGFILE.fits DIFFUSE_MASK.fits"
    #print "USE: cont_seg_all.py Ha_map.fits FLUX_LIMIT_PEAK MAX_DIST_PEAK FRAC_PEAK MIN_FLUX SEGFILE.fits DIFFUSE_MASK.fits"
    sys.exit(0)

data=sys.argv[1].split(",")
Ha_file=data[0]
e_Ha_file=data[1]
F_max=float_(sys.argv[2])
data=sys.argv[3].split(",")
target_SN=float_(data[0])
min_SN=float_(data[1])
frac=float_(sys.argv[4])
MIN_FLUX=float_(sys.argv[5])
segfile=sys.argv[6]
maskfile=sys.argv[7]

[Ha_map, hdr]=gdata(Ha_file, 0, header=True)
Ha_map[np.where(np.isfinite(Ha_map)==False)]=-1e12
[e_Ha_map, hdr]=gdata(e_Ha_file, 0, header=True)
e_Ha_map[np.where(np.isfinite(e_Ha_map)==False)]=1e20
#
# We derive the radius
#
pdl_SN=Ha_map/e_Ha_map
pdl_SN_t=np.zeros(pdl_SN.shape)
pdl_SN_t[:,:]=pdl_SN[:,:]
pdl_SN_t[np.where(pdl_SN_t==0)]=1
pdl_max_dist=target_SN/pdl_SN_t
pdl_max_dist[np.where(np.isfinite(pdl_max_dist)==False)]=0

[nx,ny]=Ha_map.shape
seg_map=np.zeros([nx,ny])
mask_map=np.ones([nx,ny])
r_max_abs=0.15*np.sqrt(nx**2.0+ny**2.0)
val_max_dist=0
max_dist=0
flux=[]
fluxl=[]
nup=0
for i in range(0, nx):
    for j in range(0, ny):
        val=Ha_map[i,j]
        if val > 0 and val < 1e30:
            flux.extend([val])
            nup=nup+1
statsa=np.array(stats(np.array(flux)))
print str(nup)+" "+str(statsa)
nul=0
for i in range(0, nx):
    for j in range(0, ny):
        val=Ha_map[i,j]
        if val > 0 and val < 1e30:
            if val < statsa[2]:
                fluxl.extend([val])
                nul=nul+1
statsa=stats(np.array(fluxl))
print str(nul)+" "+str(statsa)
map_now=Ha_map
stop=0
NR=1
JP=[]
IP=[]
F_PEAK=[]
lin_p=""
while stop == 0:
    k=0
    lower_lim=F_max
    upper_lim=1e30
    ip=-1
    jp=-1
    for i in range(0, nx):
        for j in range(0, ny):
            val=map_now[i,j]
            if val > lower_lim:
                is_maxt=is_max(i,j)
                if NR == 1:
                    is_neart=1
                else:
                    is_neart=is_near(i,j,seg_map,nx,ny)
                if val < upper_lim and val > lower_lim and is_maxt == 1:
                    lower_lim=val
                    ip=i
                    jp=j
    JP.extend([jp])
    IP.extend([ip])
    F_PEAK.extend([Ha_map[ip,jp]])
    #print "ip,jp "+str(ip)+","+str(jp)+" "+str(NR)
    #lin_p=lin_p+
    print str(NR)+","+str(ip)+","+str(jp)+","+str(NR)+","+str(val_max_dist)+","+str(max_dist)#+","
    area=1
    if ip == -1 and jp == -1:
        #print lin_p+"0"
        lin_p=""
        stop=1
    else:
        #
        # Variable max dist!
        #
        val_max_dist=pdl_max_dist[ip,jp]
        #
        # Co-variance correction
        #
        M=15.0
        area_factor=np.sqrt((15.0*3.1416*val_max_dist**2.0)/(15.0+3.1415*val_max_dist**2.0-1.0))
        max_dist=round(val_max_dist*area_factor/np.sqrt(3.1416))
        if max_dist < 1:
            max_dist=1
        if max_dist > r_max_abs:
            max_dist=int(r_max_abs)
        [seg_map,mask_map]=LOOP(ip,jp,seg_map,mask_map)
        #print str(area)
        map_now=map_now*mask_map
        NR=NR+1
for i in range(0, nx):
    for j in range(0, ny):
        val=Ha_map[i,j]
        if val >= MIN_FLUX:
            mask_map[i,j]=0
        nr=seg_map[i,j]
        if nr == 0 and val != -1e12:
            seg_map=nearest(i,j,NR,IP,JP,seg_map)
        if val == 0 or val < 0.0:
            seg_map[i,j]=0.0
if pt.exists(segfile) == False:
    wfits(segfile,seg_map,hdr)
else:
    sycall("rm "+segfile)
    wfits(segfile,seg_map,hdr)
if pt.exists(maskfile) == False:
    wfits(maskfile,mask_map,hdr)
else:
    sycall("rm "+maskfile)
    wfits(maskfile,mask_map,hdr)
