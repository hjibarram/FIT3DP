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
if len(sys.argv) < 5:
    print "USE: spec_cube_plot.py INPUT_FILE.CUBE.FITS SN_map.fits signalmap.fits noisemap.fits [WMIN,WMAX]"
    sys.exit(0)

input=sys.argv[1]
sn_map=sys.argv[2]
signal_map=sys.argv[3]
noise_map=sys.argv[4]
if len(sys.argv) == 6: 
    data=sys.argv[5].split(',')
    wmin=float_(data[0])
    wmax=float_(data[1])
else:
    wmin=0.0
    wmax=0.0

y_min=1e12
y_max=-1e12
n=0
[pdl, hdr]=gdata(input, 0, header=True)
[nz,nx,ny]=pdl.shape
crval=hdr["CRVAL3"]
cdelt=hdr["CDELT3"]
crpix=hdr["CRPIX3"]

if cdelt == 0:
    cdelt=1
wave=crval+cdelt*(np.arange(0, nz)+1-crpix)

if wmin == wmax:
    i_6300=int((5590-crval)/cdelt+crpix-1)
    i_6550=int((5680-crval)/cdelt+crpix-1)
else:
    i_6300=int((wmin-crval)/cdelt+crpix-1)
    i_6550=int((wmax-crval)/cdelt+crpix-1)

pdl_slice=pdl[i_6300:i_6550,:,:]
[nz,nx,ny]=pdl_slice.shape
delta=int((i_6550-i_6300)/10)
if np.mod(delta ,2) == 0:
    delta=delta+1
pdl_k=np.ones(delta)

pdl_res=np.zeros([nz,nx,ny])
for j in range(0, ny):
    for i in range(0, nx):
        t = pdl_slice[:,i,j]
        #print delta
        smoothed = medfilt(t,delta)
        #if np.sum(t) != 0:
        #    print smoothed,"SMOT"
        #    print t,"obs"
        pdl_res[:,i,j]=abs(t-smoothed)

map_sigma_flux=np.mean(pdl_res,axis=0)
map_mean_flux=np.mean(pdl_slice,axis=0)+np.std(pdl_slice,axis=0)

#a=np.array([1,2,3,4,5,6,5,4,3,2,1])
#print medfilt(a,9)
map_sigma_flux[np.where(np.isfinite(map_sigma_flux)==False)]=1.0
map_sigma_flux[np.where(map_sigma_flux == 0.0)]=1.0
map_mean_flux[np.where(np.isfinite(map_mean_flux)==False)]=0.0

pdl_SN=map_mean_flux/map_sigma_flux

#print np.sum(map_mean_flux)
#print np.sum(map_sigma_flux)
#print np.sum(pdl_SN)


h=pyf.PrimaryHDU().header
h["NAXIS"]=2 
h["NAXIS1"]=nx 
h["NAXIS2"]=ny

wfits(sn_map,pdl_SN,h)
wfits(signal_map,map_mean_flux,h)
wfits(noise_map,map_sigma_flux,h)
