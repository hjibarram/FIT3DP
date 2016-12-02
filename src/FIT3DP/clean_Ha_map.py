#! /usr/bin/env python
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
sys.argv=filter(None,sys.argv)
if len(sys.argv) < 2:
    print "USE: clean_Ha_map.py INPUT1.FITS INPUT2.FITS [MAX_VEL] [MIN_VEL]\n";
    sys.exit(0)
infile1=sys.argv[1]
infile2=sys.argv[2]
max_vel=1e12
if len(sys.argv) == 4:
  max_vel=float_(sys.argv[3])
min_vel=0
if len(sys.argv) == 5:
  max_vel=float_(sys.argv[3])
  min_vel=float_(sys.argv[4])
[a_in1, h ]=gdata(infile1, 0, header=True)
[a_in2, h2]=gdata(infile2, 0, header=True)
[nx,ny]=a_in1.shape
a=[]
for i in range(0, nx):
    for j in range(0, ny):
        vel=a_in1[i,j]
        mask=a_in2[i,j]
        if mask == 1:
            a.extend([vel])
pdl_a=np.array(a)
statst=my.stats(pdl_a)
for i in range(0, nx):
    for j in range(0, ny):
        vel=a_in1[i,j]
        mask=a_in2[i,j]
        if (vel-statst[2]) > 300 or vel > max_vel:
            lo=6562*(1+vel/300000.);
            new_vel=((lo/6583)-1)*300000.
            a_in1[i,j]=new_vel
            a_in2[i,j]=1
        if vel < min_vel:
            a_in2[i,j]=0
wfits(infile1, a_in1, h)
wfits(infile2, a_in2, h)