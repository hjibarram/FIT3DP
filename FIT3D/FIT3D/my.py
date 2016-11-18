#!/usr/bin/python

import sys
import numpy as np
from numpy import float_
from numpy import absolute as abs
from numpy import random as ran
import matplotlib
from scipy.signal.signaltools import convolve2d
from scipy.interpolate.interpolate import interp1d

def A_l(Rv,l):
    l=l/10000.; #Amstrongs to Microns
    x=1/l
    if x > 1.1:
        y=(x-1.82)
        ax=1+0.17699*y-0.50447*y**2-0.02427*y**3+0.72085*y**4+0.01979*y**5-0.77530*y**6+0.32999*y**7
        bx=1.41338*y+2.28305*y**2+1.07233*y**3-5.38434*y**4-0.62251*y**5+5.30260*y**6-2.09002*y**7
    else:
        ax=0.574*x**1.61
        bx=-0.527*x**1.61
    Arat=ax+bx/Rv
    return Arat

def median_filter(box,arra):
    if box == 2*int(box/2.):
        box=box+1
    val=arra
#    print val.shape,box
    for i in range(box, (len(val)+1-box)):
        tmp=np.zeros(2*box)
        for jk in range(0, 2*box):
            tmp[jk]=arra[i-box+jk]
        val[i]=np.median(tmp)
    for i in range(1, box):        
        effec_box=i
        tmp=np.zeros(2*effec_box)
        for jk in range(0, 2*effec_box):
            tmp[jk]=arra[i-effec_box+jk]
        val[i]=np.median(tmp)
    for i in range(len(val)+1-box, len(val)):        
        effec_box=len(val)-i+1
        tmp=np.zeros(2*effec_box)
        for jk in range(0, 2*effec_box):
            tmp[jk]=arra[i-effec_box+jk-1]
        val[i]=np.median(tmp)
    val[0]=val[1]
    return val

def median_box(box, arra):
    if box == 2*int(box/2.):
        box=box+1
    in_val=arra
    out_val=[]
    k=0
    for i in range(box, len(in_val)+1-box, 2*box):
        tmp=np.zeros(2*box)
        for j in range(0, 2*box):
            tmp[j]=arra[i-box+j]
        out_val.extend([np.median(tmp)])
    out_val=np.array(out_val)
    return out_val

def sycall(comand):
    from subprocess import call
    line=comand.split(" ")
    fcomand=[]
    fcomand.extend(line)
    call(fcomand)
    
def stats(data):
    out=np.zeros(7)
    out[0]=np.mean(data)
    out[1]=np.mean(data)+np.std(data)
    out[2]=np.median(data)
    out[3]=np.amin(data)
    out[4]=np.amax(data)
    out[5]=np.std(data)
    out[6]=np.mean(data)+np.std(data)
    return out
    
def mycall(comand, alias_py='python'):
    alias_py=''
    from subprocess import call
    line=comand.split(" ")
    if alias_py != '' :
        fcomand=[alias_py]
    else:
        fcomand=[]#2.7']
    fcomand.extend(line)
    linp=''
    nx=len(fcomand)
    for i in range(0, nx):
        linp=linp+fcomand[i]+" "
    print linp+" "#| PROCESSING"
    #flog.write(linp+"\n")
    #print fcomand
    call(fcomand)
    print "DONE"