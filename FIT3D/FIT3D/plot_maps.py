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

def sycallo(comand):
    import os
    out=os.popen(comand, 'r')
    line=out.readline()
    return line.replace("\n","")

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
if len(sys.argv) < 9:
    print "USE: plot_maps.pl map.fits min max bright contrast Label factor dev [MASK_MAP]";
    sys.exit(0)
mapfile=sys.argv[1]
min=float_(sys.argv[2])
max=float_(sys.argv[3])
bright=sys.argv[4]
contrast=sys.argv[5]
label=sys.argv[6]
factor=float_(sys.argv[7])
dev=sys.argv[8]
id_mask=0
if len(sys.argv) == 10:
    mask_file=sys.argv[9]
    id_mask=1
    
print dev
if dev == "null":
    import matplotlib.pyplot as plt
#    from matplotlib import gridspec
#    import matplotlib.colors as colors
#    import matplotlib.colorbar as colorbar
#    import matplotlib.cm as cmx
    import matplotlib.image as mpimg
else:
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
#    from matplotlib import gridspec
#    import matplotlib.colors as colors
#    import matplotlib.colorbar as colorbar
#    import matplotlib.cm as cmx
    import matplotlib.image as mpimg
root_dir=sycallo("echo $FIT3DP_PATH")    
table="califaCT"
reverse=0  
color_cont=0
[pdl_map, h]=gdata(mapfile, 0, header=True)
[nx,ny]=pdl_map.shape
crval2=h["CRVAL2"]
cdelt2=h["CD2_2"]
cdelt1=h["CD1_1"]
crval1=h["CRVAL1"]
#crval1=-138.839996337891
#cdelt1=4
#crval2=-210.726806640625
#cdelt2=4
if id_mask == 1:
    pdl_mask=gdata(mask_file)
    pdl_map=pdl_map*pdl_mask*factor
else:
    pdl_map=pdl_map*factor

NNx=nx-1
NNy=ny-1
rpdl_map=pdl_map
map=rpdl_map
max=np.amax(map)
min=np.amin(map[np.where(map > 0)])
#crval1=crval1-50-10
#crval2=crval2+17-10
x_max=crval1
x_min=crval1+cdelt1*nx
y_min=crval2
y_max=crval2+cdelt2*ny
for j in range(0, ny):
    for i in range(0, nx):
        if map[i,j] > 1e12:
            map[i,j]=0
        if map[i,j] == np.nan:
            map[i,j]=0
ft=open(root_dir+"/ctable", "r")
if table != "califa":
    if table == "califa2":
        g=np.zeros(254)
        r=np.zeros(254)
        b=np.zeros(254)
        l=np.zeros(254)
        for cmap in ft:
            data=cmap.split(" ")
            data=filter(None,data)
            nc=int(float_(data[0]))
            b[238-nc]=int(float_(data[1]))
            g[238-nc]=int(float_(data[2]))
            r[238-nc]=int(float_(data[3]))
            l[nc]=nc/238
        bright=1
    else:
        if table == "califaCT":
            nc=256
            g=np.zeros(nc)
            r=np.zeros(nc)
            b=np.zeros(nc)
            l=np.zeros(nc)
            for cmap in ft:
                data=cmap.split(" ")
                data=filter(None,data)
                nc=int(float_(data[0]))
                r[nc-1]=(float_(data[1]))/255.
                g[nc-1]=(float_(data[2]))/255.
                b[nc-1]=(float_(data[3]))/255.
                l[nc-1]=nc/255.
            bright=1 
            contrast=0.5
            r[0]=1.0
            g[0]=1.0
            b[0]=1.0
            nc=256
        else:
            nc=pl.shape[0]
            g=np.zeros(nc)
            r=np.zeros(nc)
            b=np.zeros(nc)
            l=np.zeros(nc)          
            for j in range(0, nc):
                l[j] = pl[j]
                r[j] = pr[j]
                g[j] = pg[j]
                b[j] = pb[j]
else:
    nc=256
    R=-1.25
    S=0.5
    gamma=0.75
    g=np.zeros(nc)
    r=np.zeros(nc)
    b=np.zeros(nc)
    l=np.zeros(nc)
    for j in range(0, nc):
        if reverse == 0:
            jj=j
        else:
            jj=255-j
        l[j] = j*(1/256.0);
        t=2*np.pi*(S/3.0+R*l[j])
        a=((l[j])**gamma)*(1-(l[j])**gamma)/2.0
        F=0.5+0.5*l[j]  
        g[jj] = F*((l[j])**gamma)-0.1486*a*cos(t)+1.7828*a*sin(t)
        r[jj] = F*((l[j])**gamma)-0.29223*a*cos(t)-0.9065*a*sin(t)
        b[jj] = F*((l[j])**gamma)+1.97294*a*cos(t)

for j in range(0, 0):
    r[j]=1
    g[j]=1
    b[j]=1

my_rgb=np.zeros([256,3])
my_rgb[:,0]=r
my_rgb[:,1]=g
my_rgb[:,2]=b
#my_rgb[:,3]=l
my_cmap = matplotlib.colors.ListedColormap(my_rgb, name='my_name')
#fig = plt.figure(figsize=(7,7))
plt.xlabel("RA (arcsec)",fontsize=14)
plt.ylabel("DEC (arcsec)",fontsize=14)
label=label.replace('Ha',r'$H_{\alpha}$')
plt.title(r''+label,fontsize=15)
plt.imshow(map,cmap=my_cmap,extent=[x_min,x_max,y_min,y_max], interpolation='none', vmin=min, vmax=max, origin='lower')
plt.colorbar()
if dev == "null":
    plt.show()
    plt.close()
else:
    plt.savefig(dev,dpi = 1000)
    plt.close()