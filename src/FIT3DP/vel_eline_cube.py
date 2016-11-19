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
if len(sys.argv) < 7:
    print "USE: vel_elines_cube.pl CUBE.fits nsearch LIMIT(% max) PREFIX_OUT WAVE_REF DEV [WMIN WMAX]"
    print "OUTPUT: PREFIX_OUT.vel_map.fits"
    print "OUTPUT: PREFIX_OUT.mask_map.fits"
    sys.exit(0)

input=sys.argv[1]
nsearch=int(sys.argv[2])
imin=float_(sys.argv[3])
outfile=sys.argv[4]
wave_ref=float_(sys.argv[5])
dev=sys.argv[6]

map_outfile=outfile+".vel_map.fits"
mask_outfile=outfile+".mask_map.fits"


if len(sys.argv) == 9:
    wmin=float_(sys.argv[7])
    wmax=float_(sys.argv[8])
else:
    wmin=0
    wmax=1e12
dmin=0

[pdl_in, hdr]=gdata(input, 0, header= True)
crval3=hdr["CRVAL3"]
cdelt3=hdr["CDELT3"]
crpix3=hdr["CRPIX3"]
[nz,nx,ny]=pdl_in.shape
map_out=np.zeros([nx,ny])
mask_out=np.zeros([nx,ny])

for ii in range(0, nx):
    for jj in range(0, ny):
# Define data
        n=0
        wave=[]
        flux=[]
        y_min=1e12
        y_max=-1e12
        t=0
        for k in range(0, nz):
            if t == 0:
                wave.extend([crval3+cdelt3*k])
                flux.extend([pdl_in[k,ii,jj]])
                t=1
            else:
                wave[n]=crval3+cdelt3*k
                flux[n]=pdl_in[k,ii,jj]
            if wave[n] > wmin and wave[n] < wmax:       
                if flux[n] > y_max:
                    y_max=flux[n]
                if flux[n] < y_min:
                    y_min=flux[n]
                t=0
                n=n+1
        # Look for peaks
        crval=wave[0]
        cdelt=wave[1]-wave[0]
        npeaks=0
        peak_y_max=[]
        peak_y_pixel=[]
        for j in range(nsearch, (n-nsearch)):
            peak=1
            for i in range(0, nsearch):
                if flux[j-i] < flux[j-i-1]:
                    peak=0
                if flux[j+i] < flux[j+i+1]:
                    peak=0
            if peak == 1:
                rrr=imin*y_max
                if flux[j] < imin*y_max:
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
                den=(fc-2*fb+fa)
                if den != 0:
                    peak_y_max.extend([c-(b-a)*((fc-fb)/den+0.5)])
                else:
                    peak_y_max.extend([0])
                npeaks=npeaks+1
        #
        # We plot the section
        #
        if y_max > y_min and y_min != 0:
            wave_peak=np.zeros(npeaks)
            #pgbegin(0,"$dev",1,1);
            #pgsfs(1.2);
            #pgscf(2);             # Set character font
            #pgslw(2);             # Set line width
            #pgsch(1.2);           # Set character height
            #pgenv($wmin,$wmax,$y_min,$y_max,0,0);
            #pglabel("Y-axis","Counts","");
            #pgline($n,\@wave,\@flux);    
            #pgsch(1.0);           # Set character height
            for k in range(0, npeaks):
                wave_peak[k]=crval+cdelt*peak_y_max[k]
                #pgsci(2);
                x=crval+cdelt*peak_y_pixel[k]
                y=0.8*y_max
                #pgpoint(1,[$x],[$y],5);
                #pgsci(4)
                x=crval+cdelt*peak_y_max[k]
                y=0.6*y_max
                #pgpoint(1,[$x],[$y],2)
            #pgsci(3);
            #pgline(2,[$wmin,$wmax],[$imin*$y_max,$imin*$y_max]);
            #pgsci(1);
            #pgsch(1.2);           # Set character height
            #pgsci(1);
            #pgclose;
            #pgend;
            if npeaks == 1:
                vel=(wave_peak[0]/wave_ref-1)*vel_light
                map_out[ii,jj]=vel
                mask_out[ii,jj]=1
            if npeaks == 2:
                vel=(wave_peak[0]/wave_ref-1)*vel_light
                map_out[ii,jj]=vel
                mask_out[ii,jj]=2
            if npeaks == 3:
                vel=(wave_peak[1]/wave_ref-1)*vel_light
                map_out[ii,jj]=vel
                mask_out[ii,jj]=2
    print str(ii)+"/"+str(nx)
hdr["NAXIS"]=2
hdr["NAXIS1"]=nx
hdr["NAXIS2"]=ny
del hdr["NAXIS3"]
hdr["COMMENT"]="vel_eline_cube result"
hdr["WMIN"]=wmin
hdr["WMAX"]=wmax
hdr["WREF"]=wave_ref
hdr["FILENAME"]=map_outfile
wfits(map_outfile, map_out, hdr)
hdr["FILENAME"]=mask_outfile
wfits(mask_outfile, mask_out, hdr)