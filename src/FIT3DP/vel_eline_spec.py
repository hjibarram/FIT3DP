#! /usr/bin/env python

import sys
import numpy as np
from numpy import float_
from numpy import absolute as abs
from numpy import random as ran
import matplotlib

vel_light=299792.458
vel=0

if len(sys.argv) <7:
    print "USE: vel_elines_cube.py spec.txt nsearch LIMIT(% max) PREFIX_OUT WAVE_REF DEV [WMIN WMAX] File_out"
    print "OUTPUT: PREFIX_OUT.vel_map.fits"
    print "OUTPUT: PREFIX_OUT.mask_map.fits"
    sys.exit(0)

input=sys.argv[1]
nsearch=int(sys.argv[2])
imin=float_(sys.argv[3])
outfile=sys.argv[4]
wave_ref=float_(sys.argv[5])
dev=sys.argv[6]
wmin=0
wmax=1e12
file_out="none"
if len(sys.argv)==9:
    wmin=float_(sys.argv[7])
    wmax=float_(sys.argv[8])
if len(sys.argv)==10:
    wmin=float_(sys.argv[7])
    wmax=float_(sys.argv[8])
    file_out=sys.argv[9]
dmin=0;


wave1=[]
spec1=[]
f=open(input,'r')
for line in f:
    if line != "\n":
        data=line.split(" ")
        data=filter(None,data)
        if data[0] != "#":
            wave1.extend([float(data[1])])
            spec1.extend([float_(data[2])])
nz=len(wave1)
f.close()

pld_wave=np.array(wave1,dtype=float_)
pdl_in=np.array(spec1,dtype=float_)

n=0
flag1=0
wave=[]
flux=[]
y_min=1e12;
y_max=-1e12;
for k in range(0, nz):
    if flag1 == 0:
        wave.extend([wave1[k]])
        flux.extend([spec1[k]])
        flag1=1
    else:
        wave[n]=wave1[k]
        flux[n]=spec1[k]
    if (wave[n]>wmin) and (wave[n] < wmax):        
        if flux[n] > y_max:
            y_max=flux[n]
        if flux[n] < y_min:
            y_min=flux[n]
        n=n+1
        flag1=0

med=np.mean(flux)
sig=np.std(flux)

crval=wave[0]
cdelt=wave[1]-wave[0]
npeaks=0
peak_y_pixel=[]
peak_y_max=[]
for j in range(nsearch, n-nsearch):
    if flux[j] > (med+2.*sig):
        peak=1
    else:
        peak=0
    for i in range(0, nsearch):
        if flux[j-i] < flux[j-i-1]:
            peak=0
        if flux[j+i] < flux[j+i+1]:
            peak=0
    if peak == 1:
        rrr=imin*y_max
        if flux[j] < (imin*y_max):
            peak=0
    if peak == 1:
        if npeaks > 0:
            delta=j-peak_y_pixel[npeaks-1]
            if delta < dmin:
                peak=0    
    if peak == 1:
        peak_y_pixel.extend([j])   
        a=j-1
        b=j
        c=j+1
        fa=-flux[a]
        fb=-flux[b]
        fc=-flux[c]
        den=fc-2*fb+fa
        if den != 0:
            peak_y_max.extend([c-(b-a)*((fc-fb)/den+0.5)])
        else:
            peak_y_max.extend([0])
        npeaks=npeaks+1
#
# We plot the section
#
wave_peak=np.zeros(npeaks)
if y_max > y_min and y_min != 0 :
    deva=dev.split("/")
    if deva[1] == "0":
        dev_plot="null"
        import matplotlib.pyplot as plt
        import matplotlib.colors as colors
        import matplotlib.cm as cmx
    else:
        dev_plot=deva[0]
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        import matplotlib.colors as colors
        import matplotlib.cm as cmx
    fig = plt.figure()
    plt.axis([wmin,wmax,y_min,y_max])
    plt.xlabel("Wavelength",fontsize=14)
    plt.ylabel("Counts",fontsize=14)
    plt.title("",fontsize=15) 
    plt.plot(wave,flux)
    for k in range(0, npeaks):
        wave_peak[k]=crval+cdelt*peak_y_max[k]
        x=crval+cdelt*peak_y_pixel[k]
        y=0.8*y_max
        poin=plt.Circle((x,y), 0.2,color='green')
        plt.gca().add_artist(poin)
        x=crval+cdelt*peak_y_max[k]
        y=0.6*y_max
        poin=plt.Circle((x,y), 0.2,color='red')
        plt.gca().add_artist(poin)
    plt.plot([wmin,wmax],[imin*y_max,imin*y_max],color='red')
    if dev_plot == "null":
        plt.show()
    else:
        plt.savefig(dev_plot)
    if npeaks == 1:
        vel=(wave_peak[0]/wave_ref-1)*vel_light
    if npeaks == 2:
        vel=(wave_peak[0]/wave_ref-1)*vel_light
    if npeaks == 3:
        vel=(wave_peak[1]/wave_ref-1)*vel_light
if file_out == "none":
    print input+", "+str(npeaks)+", "+str(vel)
else:
    print input+", "+str(npeaks)+", "+str(vel)
    f=open(file_out,'w')
    f.write(input+","+str(npeaks)+","+str(vel)+" \n")
    f.close()