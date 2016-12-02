#! /usr/bin/env python

import sys
import numpy as np
from numpy import absolute as abs
import pyfits
from pyfits import getheader as ghead
from pyfits import getdata as gdata

sys.argv=filter(None,sys.argv)

if len(sys.argv) < 3:
    print "USE: read_img_header.py FILE.FITS HEADER [OutFile]"
    sys.exit(0)

file=sys.argv[1]
list_header=sys.argv[2];
outf="none"
if len(sys.argv) == 4:
    outf=sys.argv[3]
hdr=ghead(file)
headers=list_header.split("\,")
#print hdr.keys()
values=""
f=open(outf,'w')
for i in range(0, len(headers)):
    if headers[i] in hdr:
        value=hdr[headers[i]]
    else:
        value="none"   
    values=values+" "+str(value)
    f.write(file+" "+headers[i]+" "+str(value)+" \n")
f.close()