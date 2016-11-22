#!/usr/bin/python

import sys
import numpy as np
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
    print "USE: get_slice.py INPUT_CUBE.fits PREFIX CONF_FILE"
    print "CONF_FILE: NAME START_W END_W"
    sys.exit(0)    

input_cube=sys.argv[1]
prefix=sys.argv[2]
conf_file=sys.argv[3]

name=[]
start_w=[]
end_w=[]
f=open(conf_file, "r");
for line in f:
    if line != "\n":
        data=line.split(" ")
        if data[0] != "#":
            data=filter(None,data)
            name.extend([data[0]])
            start_w.extend([float(data[1])])
            end_w.extend([float(data[2])])
f.close()
ns=len(name)
print str(ns)+" slices to cut"
print "Reading cube"
[b, hdr] = gdata(input_cube, 0 ,header=True)
crval3 = hdr["CRVAL3"]
cdelt3= hdr["CDELT3"]
crpix3= hdr["CRPIX3"]
[nx,ny,nz]=b.shape

for i in range(0, ns):
    out_file=prefix+"_"+name[i]+"_"+str(int(start_w[i]))+"_"+str(int(end_w[i]))+".fits"
    start_i=int((start_w[i]-crval3)/cdelt3+crpix3-1)
    end_i=int((end_w[i]-crval3)/cdelt3+crpix3-1)
    if (start_i > -1 ) and (end_i < nx): 
        npix=end_i-start_i+1
        a=b[start_i:end_i,:,:]
        c=np.average(a,0)
        hdr["PIX_WIDT"]=npix
        hdr["START_W"]=start_w[i]
        hdr["END_W"]=end_w[i]
        hdr["NAXIS"]=2
        del hdr["NAXIS3"]
        if pt.exists(out_file) == False:
            wfits(out_file,c,hdr)
        else:
            sycall("rm "+out_file)
            wfits(out_file,c,hdr)
        print out_file+" saved"
    else:
        print "section ("+str(int(start_w[i]))+","+str(int(end_w[i]))+") out of margings"