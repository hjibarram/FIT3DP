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
import my as my
from scipy.signal import medfilt2d
from scipy.signal import medfilt

def sycall(comand):
    from subprocess import call
    line=comand.split(" ")
    fcomand=[]
    fcomand.extend(line)
    linp=''
    nx=len(fcomand)
    for i in range(0, nx):
        linp=linp+fcomand[i]+" "
    print linp
    call(fcomand)

def wfits(name, data, hdr):
    if ptt.exists(name) == False:
        wfit(name,data,hdr)
    else:
        name1=name.replace("\ "," ")
        name1=name1.replace(" ","\ ")
        sycall("rm "+name1)
        wfit(name,data,hdr)

vel_light=299792.458
sys.argv=filter(None,sys.argv)
if len(sys.argv) < 2:
    print "USE: csv_to_map.seg.py CSVFILE column seg_file.fits map.fits"
    sys.exit(0)

infile=sys.argv[1]
nc=int(sys.argv[2])
seg_file=sys.argv[3]
map_file=sys.argv[4]
norm_file="norm_"+sys.argv[4]
area_file="area_"+sys.argv[4]

[pdl_seg, h]=gdata(seg_file, 0, header=True)
[nx,ny]=pdl_seg.shape
pdl_map=np.zeros([nx,ny])
pdl_sum=np.zeros([nx,ny])

print "Reading input files"
f=open(infile,"r")
n=1
val=[]
for line in f:
    if not "#" in line:
        data=line.replace("\n","").split(",")
        data=filter(None,data)
        id=n
        val.extend([float_(data[nc])])
        n=n+1
f.close()
val=np.array(val)
sum=np.zeros(n-1)

for i in range(0, nx):
    for j in range(0, ny):
        id=pdl_seg[i,j]
        if id > 0:
            sum[id-1]=sum[id-1]+1
        val_now=val[id-1]
        pdl_map[i,j]=val_now
        id=pdl_seg[i,j]

for i in range(0, nx):
    for j in range(0, ny):
        id=pdl_seg[i,j]
        pdl_sum[i,j]=sum[id-1]

wfits(map_file,pdl_map,h)
pdl_norm=pdl_map/pdl_sum
wfits(norm_file,pdl_norm,h)
wfits(area_file,pdl_sum,h)
print "DONE"
