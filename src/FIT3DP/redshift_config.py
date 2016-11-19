#!/usr/bin/python

import sys
import numpy as np
from numpy import float_

sys.argv=filter(None,sys.argv)

if len(sys.argv) < 3:
    print "USE: redshift_config.py input_config redshift output_config"
    sys.exit(0)

input=sys.argv[1]
z=float_(sys.argv[2])
output=sys.argv[3]

n=0
LINE=[]
f=open(input,"r")
for line in f:
    LINE.extend([line])
n=len(LINE)
f.close()

f=open(output,"w")
for i in range(0, 4):
    f.write(LINE[i])
nc=int(LINE[3])
print "NC="+str(nc)
for i in range(4, 4+nc):
    data=LINE[i].split(" ")
    data[0]=str(int(float_(data[0])*(1+z)))
    data[1]=str(int(float_(data[1])*(1+z)))
    linep=""
    for j in range(0, len(data)):
        data[j]=data[j].replace('\n','')
        linep=linep+data[j]+" "
    f.write(linep+"\n")
for i in range(4+nc,n):
    f.write(LINE[i])
f.close()

