#! /usr/bin/env python
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
from numpy import random as ran
#from aifc import data
import matplotlib
from matplotlib.path import Path
import matplotlib.patches as patches

def vertx_poly(x,y,colorn,axp):
    pols=zip(x,y)
    nl=len(x)
    codes=[]
    if nl > 2:
        codes.extend([Path.MOVETO])
        for i in range(0, nl-2):
            codes.extend([Path.LINETO])
        codes.extend([Path.CLOSEPOLY])
    else:
        print "NO sufficient vertex to plot a polygon"
        codes.extend([Path.MOVETO])
        codes.extend([Path.CLOSEPOLY])
    path = Path(pols, codes)
    patch = patches.PathPatch(path, facecolor=colorn, lw=2)
    axp.add_patch(patch)

def get_seconds():
    import time
    [year,mon,mday,hour,min,sec,wday,yday,isdst] = time.localtime()
    sec_now=yday*3600.*24.+hour*3600.+min*60.+sec
    return sec_now

def get_time():
    import time
    [year,mon,mday,hour,min,sec,wday,yday,isdst] = time.localtime()
    time="# TIME "+str(sec)+" "+str(min)+" "+str(hour)+" "+str(mday)+" "+str(mon)+" "+str(year)+" "+str(wday)+" "+str(yday)+" "+str(isdst)
    return time

def get_data_t():
    dat=[]
    dat.extend([["Ha"     ,6542.800, 6582.800, 6507.500, 6537.500, 6587.900, 6617.900]])
    dat.extend([["Hd"     ,4083.500, 4122.250, 4041.600, 4079.750, 4128.500, 4161.000]])
    dat.extend([["Hb"     ,4847.875, 4876.625, 4827.875, 4847.875, 4876.625, 4891.625]])
    dat.extend([["Mgb"    ,5160.125, 5192.625, 5142.625, 5161.375, 5191.375, 5206.375]])
    dat.extend([["Fe5270" ,5245.650, 5285.650, 5233.150, 5248.150, 5285.650, 5318.150]])
    dat.extend([["Fe5335" ,5312.125, 5352.125, 5304.625, 5315.875, 5353.375, 5363.375]])
    dat.extend([["D4000"  ,4050.000, 4250.000, 3750.000, 3950.000, 0.000   , 1.000   ]])
    dat.extend([["Hdmod"  ,4083.500, 4122.250, 4079    , 4083    , 4128.500, 4161.000]])
    dat.extend([["Hg"     ,4319.75 , 4363.50 , 4283.50 , 4319.75 , 4367.25 , 4419.75 ]])
    return dat

time1=get_seconds()


sys.argv=filter(None,sys.argv)
if len(sys.argv) < 5:
    print "USE:get_index_output_auto_ssp_elines_rnd_rss_outFIT3D.py output_auto_ssp_elines_several.RSS.cube.fits NSIM auto_ssp.rss.out DEV\n";    
    sys.exit(0)

spec_file=sys.argv[1]
NSIM=int(sys.argv[2])
file_z=sys.argv[3]
dev=sys.argv[4]
if len(sys.argv) == 5:
    printo=1
if len(sys.argv) == 6:
    printo=0
    file_out=sys.argv[5]
    f2=open(file_out, "w")

plot=0
if dev != "/null":
    if dev == "window":
        import matplotlib.pyplot as plt
        import matplotlib.colors as colors
        import matplotlib.cm as cmx
        plot=1
    else:
        matplotlib.use('Agg')
        import matplotlib.backends.backend_pdf
        import matplotlib.pyplot as plt
        import matplotlib.colors as colors
        import matplotlib.cm as cmx
        pdf = matplotlib.backends.backend_pdf.PdfPages(dev)
        plot=2

n=0
z=[]
f=open(file_z, "r")
for line in f:
    data=line.split(",")
    data=filter(None,data)
    if not "#" in data[0]:
        z.extend([float_(data[7])])    
f.close()
n=len(z)

w_min=1e12
w_max=-1e12
dat=get_data_t()
ni=len(dat)
name=[]
OL1=np.zeros(ni)
OL2=np.zeros(ni)
OLb1=np.zeros(ni)
OLb2=np.zeros(ni)
OLr1=np.zeros(ni)
OLr2=np.zeros(ni)
OLb=np.zeros(ni)
OLr=np.zeros(ni)
for i in range(0, ni):
    data=dat[i]
    name.extend([data[0]])
    OL1[i]=data[1]
    OL2[i]=data[2]
    OLb1[i]=data[3]
    OLb2[i]=data[4]
    OLr1[i]=data[5]
    OLr2[i]=data[6]
    OLb[i]=(OLb1[i]+OLb2[i])/2.
    OLr[i]=(OLr1[i]+OLr2[i])/2.
    if OLb1[i] < w_min:
        w_min=OLb1[i]
    if OLr2[i] > w_max:
        w_max=OLr2[i]

w_min=w_min-200
w_max=w_max+200

[pdl, hdr]=gdata(spec_file, 0, header=True)
[ny,nx,nz]=pdl.shape
print ny,nx,nz
crval=hdr["CRVAL1"]
cdelt=hdr["CDELT1"]
crpix=hdr["CRPIX1"]
if crval == 0:
    crval=hdr["CRVAL3"]
    cdelt=hdr["CDELT3"]
    crpix=hdr["CRPIX3"]
fwhm=2.345*3*cdelt
nsim=NSIM

w=np.zeros(nx)
for i in range(0, nx):
    w[i]=crval+cdelt*(i+1-crpix)
coun=0
L1=np.zeros(ni)
L2=np.zeros(ni)
Lb1=np.zeros(ni)
Lb2=np.zeros(ni)
Lr1=np.zeros(ni)
Lr2=np.zeros(ni)
Lb=np.zeros(ni)
Lr=np.zeros(ni)
for j in range(0, ny):
    gas=pdl[j,:,3]
    flux_res=pdl[j,:,4]
    flux_org=pdl[j,:,0]-(gas-flux_res)
    statst=my.stats(flux_org)
    if statst[0]==0 and statst[3]==statst[4]:
        flux_org=np.ones(nx)
        flux_res=100*np.ones(nx)
    res=np.array(flux_res)
    med_res=my.median_filter(int(3.0*fwhm/cdelt),res)
    pdl_res=np.array(med_res)
    # Modified by redshift
    w_min=1e12
    w_max=-1e12
    for I in range(0, ni):
        L1[I]=OL1[I]*(1+z[j])
        L2[I]=OL2[I]*(1+z[j])
        Lb1[I]=OLb1[I]*(1+z[j])
        Lb2[I]=OLb2[I]*(1+z[j])
        Lr1[I]=OLr1[I]*(1+z[j])
        Lr2[I]=OLr2[I]*(1+z[j])
        Lb[I]=(Lb1[I]+Lb2[I])/2.0
        Lr[I]=(Lr1[I]+Lr2[I])/2.0
        if Lb1[I] < w_min:
            w_min=Lb1[I]
        if Lr2[I] > w_max:
            w_max=Lr2[I]
    w_min=w_min-200
    w_max=w_max+200
    VAL_EW=np.zeros([ni,nsim])
    med_S=np.median(flux_org)
    sig_N=np.std(flux_org)
    for k in range(0, nsim):
        pdl_noise=ran.randn(nx)
        pdl_noise=pdl_noise*pdl_res
        pdl_now=flux_org+pdl_noise
        f=pdl_now
        median=np.median(f)
        min=-0.1*abs(median)
        max=3.5*median
        #
        # Plot and analysis!
        #
        if plot != 0:
            fig = plt.figure(coun)
            coun=coun+1
            axp = fig.add_subplot(111)
            plt.axis([w_min,w_max,min-0.2*abs(min),max+0.2*abs(max)])
            plt.xlabel("Wavelength",fontsize=14)
            plt.ylabel("Flux",fontsize=14)
            jet = plt.get_cmap('jet') 
            cNorm  = colors.Normalize(vmin=0, vmax=nz-1)
            scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)
            plt.bar(w, f, align="center")    
        cdelt=w[1]-w[0]
        for J in range(0, ni):
            if name[J] != "D4000":
                Sb=0;
                Sr=0
                nr=0
                nb=0
                I1_b=(Lb1[J]-w[0]-0.5*cdelt)/cdelt
                I2_b=(Lb2[J]-w[0]-0.5*cdelt)/cdelt
                for I in range(int(I1_b+1), int(I2_b)):
                    Sb=Sb+f[I]*cdelt
                    nb=nb+1
                I=int(I1_b)
                ff=I+1-I1_b
                Sb=Sb+(f[I])*ff*cdelt
                I=int(I2_b)
                ff=I2_b-I
                Sb=Sb+(f[I])*ff*cdelt
                Sb=Sb/(Lb2[J]-Lb1[J])
                I1_r=(Lr1[J]-w[0]-0.5*cdelt)/cdelt
                I2_r=(Lr2[J]-w[0]-0.5*cdelt)/cdelt
                for I in range(int(I1_r+1),int(I2_r)):
                    Sr=Sr+f[I]*cdelt
                    nr=nr+1
                I=int(I1_r)
                ff=I+1-I1_r
                Sr=Sr+(f[I])*ff*cdelt
                I=int(I2_r)
                ff=I2_r-I
                Sr=Sr+(f[I])*ff*cdelt
                Sr=Sr/(Lr2[J]-Lr1[J])
                EW=0
                K=0
                CK=[]
                wk=[]
                I1=(L1[J]-w[0]-0.5*cdelt)/cdelt
                I2=(L2[J]-w[0]-0.5*cdelt)/cdelt
                for I in range(int(I1+1), int(I2)):
                    C=Sb*((Lr[J]-w[I])/(Lr[J]-Lb[J]))+Sr*((w[I]-Lb[J])/(Lr[J]-Lb[J]))
                    EW=EW+(1-f[I]/C)*(w[I]-w[I-1])
                    CK.extend([C])
                    wk.extend([w[I]])
                    K=k+1
                I=int(I1)
                ff=I+1-I1
                C=Sb*((Lr[J]-w[I])/(Lr[J]-Lb[J]))+Sr*((w[I]-Lb[J])/(Lr[J]-Lb[J]))
                EW=EW+(1-f[I]/C)*(w[I]-w[I-1])*ff
                I=int(I2)
                ff=I2-I
                C=Sb*((Lr[J]-w[I])/(Lr[J]-Lb[J]))+Sr*((w[I]-Lb[J])/(Lr[J]-Lb[J]))
                EW=EW+(1-f[I]/C)*(w[I]-w[I-1])*ff
                EW=EW/(1+z[j])
            else:
                Sb=0
                nb=0
                I1_b=(Lb1[J]-w[0]-0.5*cdelt)/cdelt
                I2_b=(Lb2[J]-w[0]-0.5*cdelt)/cdelt
                for I in range(int(I1_b+1), int(I2_b)):
                    Sb=Sb+f[I]*cdelt
                    nb=nb+1
                I=int(I1_b)
                ff=I+1-I1_b
                Sb=Sb+(f[I])*ff*cdelt
                I=int(I2_b)
                ff=I2_b-I
                Sb=Sb+(f[I])*ff*cdelt
                Sb=Sb/(Lb2[J]-Lb1[J])
                S=0
                K=0
                I1=(L1[J]-w[0]-0.5*cdelt)/cdelt
                I2=(L2[J]-w[0]-0.5*cdelt)/cdelt
                for I in range(int(I1+1), int(I2)):
                    S=S+f[I]*cdelt
                    K=K+1
                I=int(I1)
                ff=I+1-I1
                S=S+(f[I])*ff*cdelt
                I=int(I2)
                ff=I2-I
                S=S+(f[I])*ff*cdelt
                S=S/(L2[J]-L1[J])
                if Sb != 0:
                    EW=S/Sb
                else:
                    EW=1e16
            VAL_EW[J][k]=EW
            if plot != 0:
                vertx_poly([Lb1[J],Lb1[J],Lb2[J],Lb2[J],Lb1[J]],[max*0.5,max*0.8,max*0.8,max*0.5,max*0.5],"orange",axp)
                vertx_poly([Lr1[J],Lr1[J],Lr2[J],Lr2[J],Lr1[J]],[max*0.5,max*0.8,max*0.8,max*0.5,max*0.5],"orange",axp)
                color_val = scalarMap.to_rgba(k)
                plt.plot(wk,CK,color=color_val)
                plt.plot([Lb[J]],[Sb],'ro')
                plt.plot([Lr[J]],[Sr],'ro')
                vertx_poly([L1[J],L1[J],L2[J],L2[J],L1[J]],[max*0.5,max*0.8,max*0.8,max*0.5,max*0.5],"orange",axp)
        if plot == 1:
            plt.show()
        else:
            if plot == 2:
                pdf.savefig(fig)
                plt.close(coun-1)    
    lip=""
    for J in range(0, ni):
        lip=lip+name[J]+" " 
        tmp_EW=[]
        for k in range(0, nsim):    
            tmp_EW.extend([VAL_EW[J][k]])
        EW_now=np.mean(tmp_EW)
        s_EW_now=np.std(tmp_EW)
        lip=lip+str(EW_now)+" "+str(s_EW_now)+" " 
    lip=lip+"SN  "+str(med_S)+"  "+str(sig_N)
    #print lip
    if printo == 1:
        print lip
    else:
        f2.write(lip+" \n")
if plot == 2:
    pdf.close()
if printo == 0:
    f2.close()