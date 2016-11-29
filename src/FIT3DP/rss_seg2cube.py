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

def stats(data):
    out=np.zeros(7)
    out[0]=np.mean(data)
    out[1]=np.mean(data)+np.std(data)
    out[2]=np.median(data)
    out[3]=np.amin(data)
    out[4]=np.amax(data)
    out[5]=np.std(data)
    out[6]=np.mean(data)+np.std(data)
    return out

sys.argv=filter(None,sys.argv)
if len(sys.argv) < 4:
    print "USE: rss_seg2cube.py INPUT.RSS.fits SEGMENTATION.fits OUTPUT.CUBE.FITS"
    sys.exit(0)

infile=sys.argv[1]
segfile=sys.argv[2]
outfile=sys.argv[3]
pt=outfile
pt=pt.replace('fits','pt.txt')
[a_in, h]=gdata(infile, 0, header=True)
[N,nz]=a_in.shape
crval=h["CRVAL1"]
cdelt=h["CDELT1"]
crpix=h["CRPIX1"]
[s_in,h2]=gdata(segfile,0,header=True)
[nx,ny]=s_in.shape
out_cube=np.zeros([nz,nx,ny])
h2["CRVAL3"]=crval
h2["CDELT3"]=cdelt
h2["CRPIX3"]=crpix
h2["NAXIS"]=3
h2["NAXIS"]=nz

for i in range(0, nx):
    for j in range(0, ny):
        ist=int(s_in[i,j])
        if ist > 0:
            is_out=ist-1
            out_cube[:,i,j]=a_in[is_out,:]

if ptt.exists(outfile) == False:
    wfits(outfile,out_cube,h2)
else:
    sycall("rm "+outfile)
    wfits(outfile,out_cube,h2)

print outfile+" created"