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

def sycall(comand):
    from subprocess import call
    line=comand.split(" ")
    fcomand=[]
    fcomand.extend(line)
    call(fcomand)

sys.argv=filter(None,sys.argv)
if len(sys.argv) < 3:
    print "USE: index_seg_cube.py indices.out SEGMENTATION.fits indeces.cube.FITS"
    sys.exit(0)

infile=sys.argv[1]
segfile=sys.argv[2]
outfile=sys.argv[3]

nseg=0
n_index=0
f=open(infile, "r")
for line in f:
    if "\n" in line:
        line=line.replace('\n','')
    data=line.split(" ")
    data=filter(None,data)
    if nseg == 0:
        n_index=(len(data)+1)/3
        if n_index > 0:
            INDEX_0=np.zeros([n_index,1])
            name=["" for x in range(n_index*2)]
        else:
            INDEX_0=np.zeros([1,1])
            name=["" ]
        for i in range(0, n_index):
            name[i]=data[3*i]
            name[i+n_index]="e_"+data[3*i]
    if 3*n_index == len(data):
        if nseg == 0:
            INDEX=INDEX_0
            e_INDEX=INDEX_0
        else:
            INDEX=np.concatenate((INDEX,INDEX_0),axis=1)
            e_INDEX=np.concatenate((e_INDEX,INDEX_0),axis=1)
        for i in range(0, n_index):
            INDEX[i][nseg]=data[3*i+1]
            e_INDEX[i][nseg]=data[3*i+2]
    else:
        if nseg == 0:
            INDEX=INDEX_0
            e_INDEX=INDEX_0
        else:
            INDEX=np.concatenate((INDEX,INDEX_0),axis=1)
            e_INDEX=np.concatenate((e_INDEX,INDEX_0),axis=1)
        for i in range(0, n_index):
            INDEX[i][nseg]=0
            e_INDEX[i][nseg]=0
    nseg=nseg+1
f.close()

print INDEX

[s_in, h]=gdata(segfile, 0, header=True)
[nx,ny]=s_in.shape
nz=2*n_index
pdl_cube=np.zeros([nz,nx,ny])


for i in range(0, nx):
    for j in range(0, ny):
        ist=int(s_in[i,j])
        is_out=ist-1
        if is_out >= 0:
            for k in range(n_index):
                kk=k+n_index
                val=INDEX[k][is_out]
                e_val=e_INDEX[k][is_out]
                pdl_cube[k,i,j]=val
                pdl_cube[kk,i,j]=e_val
h["NAXIS"]=3
h["NAXIS1"]=nx
h["NAXIS2"]=ny
h["NAXIS3"]=nz
h["COMMENT"]="FIT-header"
h["FILENAME"] = outfile
for k in range(0, n_index):
    index="INDEX"+str(k)
    h[index]=name[k]
    kk=k+n_index
    index="INDEX"+str(kk)
    h[index]=name[kk]

if ptt.exists(outfile) == False:
    wfits(outfile,pdl_cube,h)
else:
    sycall("rm "+outfile)
    wfits(outfile,pdl_cube,h)