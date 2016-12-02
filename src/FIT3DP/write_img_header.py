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
import my
from my import mycall

def sycall(comand):
    from subprocess import call
    line=comand.split(" ")
    fcomand=[]
    fcomand.extend(line)
    call(fcomand)

def wfits(name, data, hdr):
    if ptt.exists(name) == False:
        wfit(name,data,hdr)
    else:
        sycall("rm "+name)
        wfit(name,data,hdr)

sys.argv=filter(None,sys.argv)
if len(sys.argv) < 2:
    print "USE: write_img_header.py FILE.FITS HEADER VALUE"
    sys.exit(0)
file=sys.argv[1]
header=sys.argv[2]
value=float_(sys.argv[3])
[pdl,h]=gdata(file, 0, header=True)
h[header]=value
wfits(file,pdl,h)