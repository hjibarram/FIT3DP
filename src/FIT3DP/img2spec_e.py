#! /usr/bin/env python

import sys
import numpy as np
from numpy import absolute as abs
import pyfits
from pyfits import getdata as gdata

sys.argv=filter(None,sys.argv)

if len(sys.argv) < 2:
    print "USE: img2spec.py INPUT_FILE.FITS NY OUTPUTFILE.txt"
    sys.exit(0)

inputf=sys.argv[1]
e_input="e_"+sys.argv[1]
NY=int(sys.argv[2])
output=sys.argv[3]
y_min=1e12
y_max=-1e12
n=0
print inputf
[pdl, hdr] = gdata(inputf, 0 ,header=True)
crval = hdr["CRVAL1"]
cdelt = hdr["CDELT1"]
crpix = hdr["CRPIX1"]
[nx,ny]=pdl.shape
pdl_e = gdata(e_input, 0 )
if cdelt == 0:
    cdelt=1
wave=[]
flux=[]
e_flux=[]
f=open(output,'w')
for i in range(0, nx):
    k=i+1
    wave.extend([crval+cdelt*(i+1-crpix)])
    flux.extend([pdl[i,NY]])
    if np.isnan(flux[i]):
        flux[i]=0
    e_flux.extend([pdl_e[i,NY]])
    if np.isnan(e_flux[i]):
        e_flux[i]=1e12
    if e_flux[i] > np.sqrt(abs(flux[i])):
        e_flux[i]=np.sqrt(abs(flux[i]))
    f.write(str(k)+" "+str(wave[i])+" "+str(flux[i])+" "+str(e_flux[i])+" \n")
f.close()
