#! /usr/bin/env python
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
    print "USE: spec_extract_cube_mean.py INPUT.CUBE.fits SEGMENTATION.fits OUTPUT.RSS.FITS"
    sys.exit(0)
infile=sys.argv[1]
segfile=sys.argv[2]
outfile=sys.argv[3]
pt=outfile
pt=pt.replace('fits','pt.txt')
[a_in, h]=gdata(infile, 0, header=True)
[nz,nx,ny]=a_in.shape
crval=h["CRVAL3"]
cdelt=h["CDELT3"]
crpix=h["CRPIX3"]
[s_in,h2]=gdata(segfile,0,header=True)
[NX,NY]=s_in.shape
inverse_s=np.zeros([NX,NY])
if nx != NX or ny != NY:
    print "Dimensions does not match ("+str(nx)+","+str(ny)+") != ("+str(NX)+","+str(NY)+")"
    sys.exit(0)
DATA=stats(s_in)
ns=int(DATA[4])
a_out=np.zeros([ns,nz])
x=np.zeros(ns)
y=np.zeros(ns)
npt=np.zeros(ns)
nd=0
for i in range(0, nx):
    for j in range(0, ny):
        ist=int(s_in[i,j])
        if ist == 0:
            inverse_s[i,j]=1
        is_out=ist-1
        s=a_in[:,i,j]
        s[np.isnan(s)]=0
        a_out[is_out,:]= a_out[is_out,:]+s    
        if is_out >= 0:
            x[is_out]=x[is_out]+j
            y[is_out]=y[is_out]+i
            npt[is_out]=npt[is_out]+1
        else:
            nd=nd+1

for j in range(0, ns):
    print str(j)+" "+str(int(npt[j]))
    a_out[j,:] = a_out[j,:]/npt[j]
hd=pyf.PrimaryHDU().header
hd["CRPIX1"]=crpix
hd["CRVAL1"]=crval
hd["CDELT1"]=cdelt
hd["CRPIX2"]=1
hd["CRVAL2"]=1
hd["CDELT2"]=1
if ptt.exists(outfile) == False:
    wfits(outfile,a_out,hd)
else:
    sycall("rm "+outfile)
    wfits(outfile,a_out,hd)
if ptt.exists("diffuse.fits") == False:
    wfits("diffuse.fits",inverse_s,h2)
else:
    sycall("rm "+"diffuse.fits")
    wfits("diffuse.fits",inverse_s,h2)
size=np.sqrt(nx**2+ny**2)/(2*ns)
f=open(pt,'w')
f.write( "C "+str(size)+" "+str(size)+" 0\n")
for i in range(0, ns):
    x[i]=x[i]/npt[i]
    y[i]=y[i]/npt[i]
    f.write(str(i)+" "+str(x[i])+" "+str(y[i])+" 1\n")
f.close

print str(outfile)+" and "+str(pt)+" created"