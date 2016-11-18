#!/usr/bin/python
import sys
import numpy as np
from numpy import float_
from numpy import absolute as abs
import pyfits as pyf
from pyfits import getheader as ghead
from pyfits import getdata as gdata
from pyfits import writeto as wfits
import os.path as ptt
import my
import matplotlib

sys.argv=filter(None,sys.argv)
if len(sys.argv) <3:
    print "USE: plot_output_auto_ssp_elines_several_Av_log_rss.py INPUT_FILE.FITS DEVICE [MIN MAX] [WMIN WMAX] [REF_LINES CUT]"
    sys.exit(0)
input=sys.argv[1]
dev=sys.argv[2]
deft=0
if len(sys.argv) ==5:
    min=float_(sys.argv[3])
    max=float_(sys.argv[4])
    deft=1
y_min=1e12
y_max=-1e12
n=0
[pdl, hdr]=gdata(input, 0, header= True);
[ny,nx,nz]= pdl.shape
crval=hdr["CRVAL1"]
cdelt=hdr["CDELT1"]
crpix=hdr["CRPIX1"]
if crval == 0:
    crval=hdr["CRVAL3"]
    cdelt=hdr["CDELT3"]
    crpix=hdr["CRPIX3"]
if cdelt == 0:
    cdelt=1
if dev == "null":
    import matplotlib.pyplot as plt
    import matplotlib.colors as colors
    import matplotlib.cm as cmx
else:
    matplotlib.use('Agg')
    import matplotlib.backends.backend_pdf
    import matplotlib.pyplot as plt
    import matplotlib.colors as colors
    import matplotlib.cm as cmx
    pdf = matplotlib.backends.backend_pdf.PdfPages(dev)    
wave=np.zeros([nx,nz])
flux=np.zeros([nx,nz])
for NY in range(0, ny):
    name=[]
    for k in range(0, nz):
        key="NAME"+str(k)
        name.extend([hdr[key]])
        for i in range(0, nx):
            wave[i][k]=crval+cdelt*(i+1-crpix)
            flux[i][k]=pdl[NY,i,k]
            if flux[i][k] > y_max:
                y_max=flux[i][k]
            if flux[i][k] < y_min:
                y_min=flux[i][k]
    wmin=wave[0][0]
    wmax=wave[nx-1][0]
    print str(wmin)+" "+str(wmax)+" "+str(NY)
    if len(sys.argv) == 5:
        min=float_(sys.argv[3])
        max=float_(sys.argv[4])
        deft=1
    if len(sys.argv) == 7:
        min=float_(sys.argv[3])
        max=float_(sys.argv[4])
        wmin=float_(sys.argv[5])
        wmax=float_(sys.argv[6])
        deft=1
    ref_def=0
    if len(sys.argv) == 9:
        min=float_(sys.argv[3])
        max=float_(sys.argv[4])
        wmin=float_(sys.argv[5])
        wmax=float_(sys.argv[6])
        ref_lines=sys.argv[7]
        cut_lines=sys.argv[8]
        ref_def=1
        deft=1
    flux_tmp=[]
    kk=0
    wave_now=np.zeros(nx)
    flux_now=np.zeros(nx)
    for i in range(0, nx):
        wave_now[i]=wave[i][0]
        flux_now[i]=flux[i][0]
        if wave_now[i] > wmin and wave_now[i] < wmax:
            flux_tmp.extend([flux_now[i]])
            kk=kk+1
    if deft == 1 and min != max:
        y_min=min
        y_max=max
    else:
        med=np.mean(flux_tmp)
        sigma=np.std(flux_tmp)
        y_min=(-0.25)*med
        y_max=2*med+sigma    
    if y_min == y_max:
        y_min=-1
        y_max=10
    fig = plt.figure(NY)
    plt.axis([wmin,wmax,y_min,y_max])
    plt.xlabel("Wavelength",fontsize=14)
    plt.ylabel("Flux",fontsize=14)
    jet = plt.get_cmap('jet') 
    cNorm  = colors.Normalize(vmin=0, vmax=nz-1)
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)
    for k in range(0, nz):
        for i in range(0, nx):
            wave_now[i]=wave[i][k]
            flux_now[i]=flux[i][k]
        color_val = scalarMap.to_rgba(k)
        plt.plot(wave_now,flux_now,color=color_val,label=name[k],lw=0.5)
        #plt.text(wmin+0.05*(wmax-wmin),y_max-0.05*(k+1)*(y_max-y_min),name[k])
        if k == 0:
            NNY=NY+1
            plt.text(wmin+0.85*(wmax-wmin),y_max-0.05*(k+1)*(y_max-y_min),"# "+str(NNY))
    plt.legend(loc=2)      
    if dev == "null":
        plt.show()
    else:
        pdf.savefig(fig)
        plt.close(NY)
if dev != "null":
    pdf.close()

sec=[]
nn=0
for k in range(0, nz):
    for i in range(0, nx):
        if wave[i][k] > wmin and wave[i][k] < wmax:
            sec.extend([flux[i][k]])
            nn=nn+1
median=np.median(sec)
mean=np.mean(sec)
sig=np.std(sec)
flux=mean*n*(wave[1]-wave[0])
print str(median)+" +- "+str(sig)+" (F="+str(flux)+")"
