#!/usr/bin/python
import sys
import numpy as np
from numpy import float_
from numpy import absolute as abs
from numpy import random as ran
import cosmolopy.distance as cd
import time
import my_auto_ssp_elines_rnd as ssp
import pyfits as pyf
from pyfits import getheader as ghead
from pyfits import getdata as gdata
from pyfits import writeto as wfit
from scipy.interpolate.interpolate import interp1d
import my as my
import os.path as ptt
import matplotlib
matplotlib.use('Agg')

def sycall(comand):
    import os
    os.system(comand)

def wfits(name, data, hdr):
    if ptt.exists(name) == False:
        wfit(name,data,hdr)
    else:
        sycall("rm "+name)
        wfit(name,data,hdr)
        
sys.argv=filter(None,sys.argv)
if len(sys.argv) < 4:
    print "USE: pack_results_name.py NAME pack_CS.csv PREFIX"
    sys.exit(0)
speed_of_light=299792.458
name=sys.argv[1]
pack=sys.argv[2]
prefix=sys.argv[3]
npi=0
LINE=[]
f=open(pack,"r")
for line in f:
    if not "#" in line:
        line=line.replace('NAME',name)
        LINE.extend([line])
        npi=npi+1
f.close()

for i in range(0, npi):
    line=LINE[i]
    data=line.split(",")
    if i == 0:
        [pdl_img, hdr]=gdata(data[1], 0, header=True)
        [nx,ny]=pdl_img.shape
        pdl_cube=np.zeros([npi,nx,ny])
        hdr["NAXIS"]=3
        hdr["NAXIS1"]=nx
        hdr["NAXIS2"]=ny
        hdr["NAXIS3"]=npi
        hdr["COMMENT"]="FITs header"
        head="OBJECT"
        hdr[head]=name
    else:
        pdl_img=gdata(data[1])
    pdl_cube[i,:,:] = pdl_img
    head="ID_"+str(i)
    hdr[head]=data[0]
    head="FILE_"+str(i)
    hdr[head]=data[1]
    head="DESC_"+str(i)
    hdr[head]=data[2]
    head="TYPE_"+str(i)
    hdr[head]=data[3]
    head="UNIT_"+str(i)
    hdr[head]=data[4].replace('\n','')
output=name+"."+prefix+".cube.fits"
wfits(output,pdl_cube,hdr)
call="gzip -f "+output
sycall(call)