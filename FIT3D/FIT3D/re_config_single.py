#!/usr/bin/python

import sys
import numpy as np
from numpy import float_ as flt
import time
import os.path as pt
import pyfits
from pyfits import getheader as ghead
from pyfits import getdata as gdata
from pyfits import writeto as wdata

def sycallo(comand):
    import os
    out=os.popen(comand, 'r')
    line=out.readline()
    return line.replace("\n","")

sys.argv=filter(None,sys.argv)        
if len(sys.argv) < 2:
    print "USE: re_config.py CONF DIR_CONF"
    sys.exit(0)
conf1=sys.argv[1]
dir_cnf=sys.argv[2]

root_dir=sycallo("echo $FIT3DP_PATH")
#print root_dir
root_dirn=root_dir.replace(" ","\ ")
f1=open(dir_cnf+"/"+conf1,"r")
f1a=open(dir_cnf+"/"+conf1+"_n","w")
#print  dir_cnf+"/"+conf1+"_n"
cont=0
for line in f1:
    cont=cont+1
    if cont == 4:
        val=flt(line.replace("\n",""))
    if cont >= 5 and cont <= 5+val-1:
        data=line.replace("\n","").split(" ")
        data=filter(None,data)
        #print line.replace("\n","")
        dir1=root_dir+"/"+data[3]
        dir2=root_dir+"/"+data[5]
        f1a.write(data[0]+" "+data[1]+" "+data[2]+" "+dir1+" "+data[4]+" "+dir2+" "+data[6]+" "+data[7]+"\n")
    else:
        f1a.write(line)
f1a.close()