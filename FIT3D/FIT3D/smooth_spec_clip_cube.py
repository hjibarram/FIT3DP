#!/usr/bin/python
import sys
import numpy as np
from numpy import float_
import pyfits as py
from pyfits import getheader as ghead
from pyfits import getdata as gdata
from pyfits import writeto as wfits
import os.path as pt
from scipy.ndimage.filters import median_filter as med2df
import my

def sycall(comand):
    from subprocess import call
    line=comand.split(" ")
    fcomand=[]
    fcomand.extend(line)
    call(fcomand)
    
sys.argv=filter(None,sys.argv)
if len(sys.argv) < 7:
    print "USE: smooth_spec_clip_cube.py intput_spec.cube.fits output_spec.cube.fits WIDTH NSIGMA NXMIN NXMAX"
    sys.exit(0)
input=sys.argv[1]
outfile=sys.argv[2]
box=int(sys.argv[3])
nsigma=float_(sys.argv[4])
nxmin=int(sys.argv[5])
nxmax=int(sys.argv[6])

[pdl_in, hdr]=gdata(input, 0, header=True)
[nz,nx,ny]=pdl_in.shape

for i in range(0, nx):
    pdl_now=pdl_in[:,i,:]
    pdl_cut=pdl_now[nxmin:nxmax,:]
    statst=my.stats(pdl_cut)
    med=statst[2]
    sig=statst[1]
    min=statst[3]
    max=med+nsigma*sig
    pdl_now=np.clip(pdl_now,min,max)
    pdl_smooth=med2df(pdl_now,size=(box,1), mode='reflect')
    pdl_in[:,i,:]=pdl_smooth
if pt.exists(outfile) == False:
    wfits(outfile,pdl_in,hdr)
else:
    sycall("rm "+outfile)
    wfits(outfile,pdl_in,hdr)