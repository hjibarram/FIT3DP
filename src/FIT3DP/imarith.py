#! /usr/bin/env python

import sys
import numpy as np
from numpy import float_
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
if len(sys.argv) < 4:
    print "USE: imarith.pl INPUT1.FITS OPERATOR(+,-,/,*) INPUT2.FITS OUTPUT.FITS"
    sys.exit(0)
infile1=sys.argv[1]
operator=sys.argv[2]
infile2=sys.argv[3]
outfile=sys.argv[4]
print infile1
[a_in1, h] = gdata(infile1, 0, header=True)
if not "fit" in infile2:
    a_in2=infile2
else:
    [a_in2, h2] = gdata(infile2, 0, header=True)
if operator == "+":
    a_in1=a_in1+a_in2
if operator == "-":
    a_in1=a_in1-a_in2
if operator == "/":
    a_in1=a_in1/a_in2
if operator == "*":
    a_in1=a_in1*a_in2
if pt.exists(outfile) == False:
    wfits(outfile,a_in1,h)
else:
    sycall("rm "+outfile)
    wfits(outfile,a_in1,h)