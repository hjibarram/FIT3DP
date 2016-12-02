#! /usr/bin/env python
import matplotlib
import sys
import numpy as np
from pyfits import getheader as ghead
from pyfits import getheader as gdata
from pyfits import writeto as wfit
import pyfits as pyf
from scipy.interpolate.interpolate import interp1d
import os.path as pt

def wfits_ext(name, hlist):
  if pt.exists(name) == False:
     hlist.writeto(name)
  else:
     sycall("rm "+name)
     hlist.writeto(name)

def wfits(name, data, hdr):
  if pt.exits(name) == False:
     wfit(name,data,hdr)
  else:
     sycall("rm "+name)
     wfit(name,data,hdr)

def sycall(comand):
  import os
  linp=comand
  os.system(comand)

sys.argv=filter(None,sys.argv)
name=sys.argv[1]

hdu_list=pyf.open(name)
flux_list= hdu_list['FLUX']
ivar_list= hdu_list['IVAR']
wave_list= hdu_list['WAVE']
mask_list= hdu_list['MASK']
flux=flux_list.data
wave=wave_list.data
nw=len(wave)
ns=flux.shape[1]
dl=wave[1]-wave[0]
h=flux_list.header
h1_l=pyf.PrimaryHDU(flux)
h1=h1_l.header
#print h[8]
#print h[23]
#print "qqqqqqqqqqqqqqqqqqqqqqqqqqqqqqq"
print h
#h['NAXIS1']=nw
#h['NAXIS2']=ns
#h1[7]=h[8]
h1['CRVAL3']=h['CRVAL3']
h1['CDELT3']=h['CD3_3']
h1['CRPIX3']=h['CRPIX3']
h1['CRVAL1']=h['CRVAL1']
#h1['CD2_1']=h['CD2_1']
h1['CD1_1']=h['CD1_1']
h1['CRPIX1']=h['CRPIX1']
h1['CTYPE1']=h['CTYPE1']
h1['CRVAL2']=h['CRVAL2']
h1['CD2_2']=h['CD2_2']
h1['CRPIX2']=h['CRPIX2']
h1['CTYPE2']=h['CTYPE2']
#h['EXTEND']=True
#h['SIMPLE']=True
hlist=pyf.HDUList([h1_l,ivar_list,mask_list,mask_list])
wfits_ext(name.replace('-LINCUBE.fits.gz','.cube.fits.gz'),hlist)
