#! /usr/bin/env python
import sys
import numpy as np
from numpy import float_
import pyfits as py
from pyfits import getheader as ghead
from pyfits import getdata as gdata
from pyfits import writeto as wfits
import os.path as pt

def sycall(comand):
    from subprocess import call
    line=comand.split(" ")
    fcomand=[]
    fcomand.extend(line)
    call(fcomand)
    
sys.argv=filter(None,sys.argv)
if len(sys.argv) < 3:
    print "USE: create_mask_map.pl map.fits cut mask.fits"
    sys.exit(0)
input_file=sys.argv[1]
cut=float_(sys.argv[2])
output_file=sys.argv[3]
[pdl_input, h]=gdata(input_file, 0, header=True)
[nx,ny]=pdl_input.shape
pdl_output=np.zeros([nx,ny])
for j in range(0, ny):
    for i in range(0, nx):
        val=pdl_input[i,j]
        if val > cut:
            pdl_output[i,j]=1
if pt.exists(output_file) == False:
    wfits(output_file,pdl_output,h)
else:
    sycall("rm "+output_file)
    wfits(output_file,pdl_output,h)