#! /usr/bin/env python

import sys
import numpy as np
from numpy import float_
import pyfits as py
from pyfits import getheader as ghead
from pyfits import getdata as gdata
from pyfits import writeto as wfits
import os.path as pt
from scipy.ndimage.filters import median_filter as med2df

def sycall(comand):
    from subprocess import call
    line=comand.split(" ")
    fcomand=[]
    fcomand.extend(line)
    call(fcomand)
    
sys.argv=filter(None,sys.argv)

if len(sys.argv) < 4:
    print "USE: med2df.py INPUT.FITS OUTPUT.fits X_WIDTH Y_WIDTH"
    sys.exit(0)

infile=sys.argv[1]
outfile=sys.argv[2]
dx=int(sys.argv[3])
dy=int(sys.argv[4])
[pdl, hdr] = gdata(infile, 0 ,header=True)
mpdl=med2df(pdl,size=(dx,dy), mode='reflect')
if pt.exists(outfile) == False:
    wfits(outfile,mpdl,hdr)
else:
    sycall("rm "+outfile)
    wfits(outfile,mpdl,hdr)
