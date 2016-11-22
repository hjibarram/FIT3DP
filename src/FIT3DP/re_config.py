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
if len(sys.argv) < 5:
    print "USE: re_config.py CONF1 CONF2 CONF3 DIR_CONF"
    sys.exit(0)
conf1=sys.argv[1]
conf2=sys.argv[2]
conf3=sys.argv[3]
dir_cnf=sys.argv[4]

root_dir=sycallo("echo $FIT3DP_PATH")+"/.."
print root_dir
root_dirn=root_dir.replace(" ","\ ")
f1=open(dir_cnf+"/"+conf1,"r")
f1a=open(dir_cnf+"/"+conf1+"_n","w")
f2=open(dir_cnf+"/"+conf2,"r")
f2a=open(dir_cnf+"/"+conf2+"_n","w")
f3=open(dir_cnf+"/"+conf3,"r")
f3a=open(dir_cnf+"/"+conf3+"_n","w")
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

cont=0
for line in f2:
    cont=cont+1
    if cont == 4:
        val=flt(line.replace("\n",""))
    if cont >= 5 and cont <= 5+val-1:
        data=line.replace("\n","").split(" ")
        data=filter(None,data)
        #print line.replace("\n","")
        dir1=root_dir+"/"+data[3]
        dir2=root_dir+"/"+data[5]
        f2a.write(data[0]+" "+data[1]+" "+data[2]+" "+dir1+" "+data[4]+" "+dir2+" "+data[6]+" "+data[7]+"\n")
    else:
        f2a.write(line)
f2a.close()

cont=0
for line in f3:
    cont=cont+1
    if cont == 4:
        val=flt(line.replace("\n",""))
    if cont >= 5 and cont <= 5+val-1:
        data=line.replace("\n","").split(" ")
        data=filter(None,data)
        #print line.replace("\n","")
        dir1=root_dir+"/"+data[3]
        dir2=root_dir+"/"+data[5]
        f3a.write(data[0]+" "+data[1]+" "+data[2]+" "+dir1+" "+data[4]+" "+dir2+" "+data[6]+" "+data[7]+"\n")
    else:
        f3a.write(line)
f3a.close()