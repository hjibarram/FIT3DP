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

sys.argv=filter(None,sys.argv)
if len(sys.argv) < 2:
    print "USE: get_SN_rss.py INPUT_FILE.RSS.FITS out_values.csv [WMIN,WMAX]"
    sys.exit(0)

input=sys.argv[1]
output=sys.argv[2]
if len(sys.argv) == 4: 
    data=sys.argv[3].split(',')
    wmin=float_(data[0])
    wmax=float_(data[1])
else:
    wmin=0.0
    wmax=0.0

y_min=1e12
y_max=-1e12
n=0
[pdl, hdr]=gdata(input, 0, header=True)
[ny,nx]=pdl.shape
crval=hdr["CRVAL1"]
cdelt=hdr["CDELT1"]
crpix=hdr["CRPIX1"]

if cdelt == 0:
    cdelt=1
wave=crval+cdelt*(np.arange(0, ny)+1-crpix)

if wmin == wmax:
    i_6300=int((5590-crval)/cdelt+crpix-1)
    i_6550=int((5680-crval)/cdelt+crpix-1)
else:
    i_6300=int((wmin-crval)/cdelt+crpix-1)
    i_6550=int((wmax-crval)/cdelt+crpix-1)
print "# "+str(i_6300)+":"+str(i_6550)+" "+str(wmin)+","+str(wmax)
pdl_slice=pdl[:,i_6300:i_6550]
[ny,nx]=pdl_slice.shape
delta=int((i_6550-i_6300)/10)
pdl_k=np.ones(delta)

pdl_res=np.zeros([ny,nx])
for j in range(0, ny):
    t = pdl_slice[j,:]
    smoothed = medfilt(t,delta)
    pdl_res[j,:]=abs(t-smoothed)

map_sigma_flux=np.mean(pdl_res,axis=1)
map_mean_flux=np.mean(pdl_slice,axis=1)+np.std(pdl_slice,axis=1)

map_sigma_flux[np.where(np.isfinite(map_sigma_flux)==False)]=1.0
map_sigma_flux[np.where(map_sigma_flux == 0.0)]=1.0
map_mean_flux[np.where(np.isfinite(map_mean_flux)==False)]=0.0

pdl_SN=map_mean_flux/map_sigma_flux

f=open(output,"w")
f.write("# (1) id\n")
f.write("# (2) S/N\n")
f.write("# (3) Signal\n")
f.write("# (4) Noise\n")
for j in range(0, ny):
    SN=pdl_SN[j]
    signal=map_mean_flux[j]
    noise=map_sigma_flux[j]
    f.write(str(j)+","+str(SN)+","+str(signal)+","+str(noise)+"\n")
f.close()
