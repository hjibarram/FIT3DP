#! /usr/bin/env python
import matplotlib
import sys
import numpy as np
from pyfits import getheader as ghead
from pyfits import getdata as gdata
from pyfits import writeto as wfit
import pyfits as pyf
from scipy.interpolate.interpolate import interp1d
import os.path as pt

def wfits(name, data, hdr):
    if pt.exists(name) == False:
        wfit(name,data,hdr)
    else:
        name1=name.replace("\ "," ")
        name1=name1.replace(" ","\ ")
        sycall("rm "+name1)
        wfit(name,data,hdr)

def sycall(comand):
    import os
    linp=comand
    os.system(comand)

def sycallo(comand):
    import os
    out=os.popen(comand, 'r')
    line=out.readline()
    return line.replace("\n","")

sys.argv=filter(None,sys.argv) 

file=sys.argv[1]
out_name=sys.argv[2]
name_a=[]
f=open(file,"r")
for line in f:
    if not "#" in line:
        name_a.extend([line.replace("\n","").replace(" ","")])

wave_min=3800.0
wave_max=9100.0
dl=1.25
wave=np.arange(wave_min,wave_max,dl)
name_a=np.array(name_a)
nw=len(wave)
ns=len(name_a)
flux_a=np.zeros([ns,nw])
flux_e=np.zeros([ns,nw])
for i in range(0, ns):
    hdu_list = pyf.open(name_a[i])
    table_hdu = hdu_list[1]
    table_hdu1 = hdu_list[2]
    table_data = table_hdu.data
    table_data1 = table_hdu1.data
    flux=table_data.field('flux')
    fluxe=table_data.field('ivar')
    logl=table_data.field('loglam')
    red=table_data1.field('Z')[0]
    fluxl = interp1d(10**logl, flux/10.0,kind='linear',bounds_error=False)(wave)
    fluxl_e = interp1d(10**logl, fluxe,kind='linear',bounds_error=False)(wave)
    flux_a[i,:] = fluxl
    flux_e[i,:] = fluxl_e
h=pyf.PrimaryHDU().header
h["NAXIS"]=2 
h["NAXIS1"]=nw 
h["NAXIS2"]=ns
h["BUNIT"]="Data Value"
h["CRVAL1"]=wave[0]
h["CDELT1"]=wave[1]-wave[0];
h["CRPIX1"]=1
wfits("CS."+out_name+".RSS.fits", flux_a, h)
wfits("e_CS."+out_name+".RSS.fits", flux_e, h)    