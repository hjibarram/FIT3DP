#!/usr/bin/python
import sys
import numpy as np
from numpy import float_
from numpy import absolute as abs
import pyfits as pyf
from pyfits import getheader as ghead
from pyfits import getdata as gdata
from pyfits import writeto as wfits
import os.path as ptt
import my

sys.argv=filter(None,sys.argv)
if len(sys.argv) < 2:
    print "USE: get_CS_slice.pl name"
    sys.exit(0)
name=sys.argv[1]
segfile="cont_seg."+name+".fits.gz"
mapHa=name+".V.fits.gz"
output="CS."+name+".slice"

[seg, h1]=gdata(segfile, 0, header=True)
[Ha, h2]=gdata(mapHa, 0, header=True)
[nx,ny]=seg.shape
n=0
ns=int(np.amax(seg))
flux=np.zeros(ns)
X=np.zeros(ns)
Y=np.zeros(ns)
A=np.zeros(ns)
for i in range(0, nx):
    for j in range(0, ny):
        ID=int(seg[i,j])
        if ID > 0:
            if n < ID:
                n=ID
            val=Ha[i,j]
            if val == 0:
                val=1
            flux[ID-1]=flux[ID-1]+val
            X[ID-1]=X[ID-1]+i*val
            Y[ID-1]=Y[ID-1]+j*val
            A[ID-1]=A[ID-1]+1

f=open(output, "w");
f.write("# (1) ID\n")
f.write("# (2) X\n")
f.write("# (3) Y\n")
f.write("# (4) Flux\n")
f.write("# (5) Area\n")
for i in range(0, n):
    if flux[i] > 0:
        x=X[i]/flux[i]
        y=Y[i]/flux[i]
    else:
        x=X[i]/flux[i]
        y=Y[i]/flux[i]
    f.write("$i  $x  $y   $flux[$i] $A[$i]\n")
f.close()