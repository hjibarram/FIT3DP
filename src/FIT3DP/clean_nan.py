#! /usr/bin/env python

import sys
import numpy as np
import pyfits
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

if len(sys.argv) < 2:
    print "USE: clean_nan.pl INPUT.FITS BADVAL"
    sys.exit(0)
    
infile1=sys.argv[1]
val=int(sys.argv[2])
[a_in1,head]=gdata(infile1,0,header=True)
a_in1[np.isnan(a_in1)]=val
if pt.exists(infile1) == False:
    wfits(infile1,c,hdr)
else:
    sycall("rm "+infile1)
    wfits(infile1,a_in1,head)