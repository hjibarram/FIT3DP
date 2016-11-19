#!/usr/bin/python

import sys
import numpy as np
from numpy import float_
from numpy import absolute as abs
import pyfits
from pyfits import getdata as gdata
from pyfits import writeto as wfits
import os.path as pt
import matplotlib

def sycall(comand):
    from subprocess import call
    line=comand.split(" ")
    fcomand=[]
    fcomand.extend(line)
    call(fcomand)

sys.argv=filter(None,sys.argv)

if len(sys.argv) < 6:
    print "USE: radial_sum_cube.pl CUBE.fits Delta_R X_C Y_C OUTPUT.RSS.fits [PLOT]";
    print "It will use images coordinates: 0,nx 0,ny"
    sys.exit(0)
    
input=sys.argv[1]
Dr=float(sys.argv[2])
x_c=int(sys.argv[3])
y_c=int(sys.argv[4])
output=sys.argv[5]
e_output="e_"+output
weight_output="weight."+output
plot=0
if len(sys.argv) == 7:
    plot=int(sys.argv[6])

[pdl_cube,hd]=gdata(input, 0, header=True)
[NX,NY,NZ]=pdl_cube.shape
y_c=NY-y_c

crval=hd["CRVAL3"]
cdelt=hd["CDELT3"]
crpix=hd["CRPIX3"]

pdl_mask=np.ones([NX,NY,NZ])                                                                                                                             
extend=hd["EXTEND"]
print "EXTEND = "+str(extend)
if (extend == 1) or (extend == "T"):
    mask_ext=input
    [b_mask,hd_m]=gdata(mask_ext, 3, header=True)
    pdl_mask=pdl_mask#-b_mask checar
    e_ext=input
    [pdl_e,hd_e]=gdata(e_ext, 1, header=True)
    w_ext=input
    [pdl_w,hd_w]=gdata(w_ext, 2, header=True)
    
pdl_cube=pdl_cube*pdl_mask
nx_med=int(NX/2)
start_i=nx_med-50
end_i=nx_med+50
a=pdl_cube[start_i:end_i,:,:]
d=np.average(a,0)
a_mask=pdl_mask[start_i:end_i,:,:]
d_mask=np.average(a_mask,0)
pdl_weight=d/d_mask
wave=[]
for i in range(0, NX):
    wave.extend([crval+cdelt*(i+1-crpix)])
wave=np.array(wave,dtype=float_)

wmin=wave[0]
wmax=wave[NX-1]

shaper="R"
size=1

nx=NX
ny=NZ*NY
pdl_in=np.zeros([nx,ny])
pdl_in_e=np.zeros([nx,ny])
pdl_in_mask=np.zeros([nx,ny])
pdl_in_weight=np.zeros(ny)
r_max=0
n=0
x_min=10e10
x_max=-10e10
y_min=10e10
y_max=-10e10
id=[]
x=[]
y=[]
r=[]
nk=[]
c=[]
s=[]
S=[]
for j in range(0, NY):
    for i in range(0, NZ):
        px=int(i+(NY-1-j)*NZ)
        pdl_in[:,px] = pdl_cube[:,j,i]
        pdl_in_e[:,px] = pdl_e[:,j,i]
        pdl_in_mask[:,px] = pdl_mask[:,j,i] 
        pdl_in_weight[px] = pdl_weight[j,i]
        id.extend([i+j*nx+1])
        x.extend([i])
        y.extend([j])
        r.extend([np.sqrt((x[n]-x_c)**2+(y[n]-y_c)**2)]) 
        nk=r[n]/Dr
        c.extend([2+int(nk/2)])
        s.extend([3+int(nk/2)])
        S.extend([2.5-nk/15])
        if r_max < r[n]:
            r_max=r[n]
        if x_min > x[n]:
            x_min=x[n]
        if y_min > y[n]:
            y_min=y[n]
        if x_max < x[n]:
            x_max=x[n]
        if y_max < y[n]:
            y_max=y[n]
        n=n+1
nr=int(r_max/Dr)+1
pdl_out=np.zeros([nx,nr])
pdl_out_e=np.zeros([nx,nr])
pdl_out_mask=np.zeros([nx,nr])
NN=[]
for i in range(0, nr):
    r_min=Dr*i
    r_max=Dr*(i+1)
    nsum=0
    sum_all=0
    t=pdl_out[:,i]
    t_e=pdl_out_e[:,i]
    t_m=pdl_out_mask[:,i]
    spaxels=[]
    nspaxels=0
    flux_max=-1e20
    flux_min=1e20
    nx1=int(0.25*nx)
    nx2=int(0.6*nx)
    slice=pdl_in[nx1:nx2,:]
    image = np.sum(slice,0)
    flux=[]
    for ns in range(0, ny):
        flux.extend([image[ns]])
        if np.isnan(flux[ns]):
            flux[ns]=0
        if flux_min > flux[ns]:
            flux_min=flux[ns]
        if flux_max < flux[ns]:
            flux_max=flux[ns]
    if np.isinf(flux_min):
        flux_min=0
    for j in range(0, ny):
        if r[j] < r_max:
            spaxels.extend([j])
            nspaxels=nspaxels+1
            slice=pdl_in[:,j]
            slice_e=pdl_in_e[:,j]
            slice_mask=pdl_in_mask[:,j]
            weight=pdl_in_weight[j]
            if flux[j] > 0:
                pdl_out[:,i] = t+slice
                pdl_out_e[:,i] = t_e+slice_e
                pdl_out_mask[:,i] = t_m+slice_mask*weight
                nsum=nsum+weight
    NN.extend([nsum])
    if NN[i] > 0 :
        pdl_out[:,i] = (t/t_m)*NN[i]
        pdl_out_e[:,i] = (t_e/t_m)*np.sqrt(abs(NN[i]))
        pdl_out_mask[:,i] = t_m/NN[i]
    if plot > 0:
        if plot == 1:
            dev="null"
            import matplotlib.pyplot as plt
            from matplotlib import gridspec
            import matplotlib.colors as colors
            import matplotlib.colorbar as colorbar
            import matplotlib.cm as cmx
        else:
            matplotlib.use('Agg')
            import matplotlib.pyplot as plt
            from matplotlib import gridspec
            import matplotlib.colors as colors
            import matplotlib.colorbar as colorbar
            import matplotlib.cm as cmx
            dev="radial_sum_cube_"+str(i)+".pdf"
        fig = plt.figure(figsize=(24,7))
        gs = gridspec.GridSpec(1, 3, width_ratios=[1.3,0.14, 3])
        ax0 = plt.subplot(gs[0])
        jet = plt.get_cmap('jet')
        n_max=100
        n_min=0 
        cNorm  = colors.Normalize(vmin=n_min, vmax=n_max)
        scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)
        plt.axis([x_min,x_max,y_min,y_max])
        plt.xlabel("RA (arcsec)",fontsize=14)
        plt.ylabel("DEC (arcsec)",fontsize=14)
        plt.title("",fontsize=15)
        for ii in range(0, ns):
            color=n_min+np.ceil(((flux[ii]-flux_min)/(flux_max-flux_min))*(n_max-n_min))
            if np.isnan(color):
                color=0
            color_val = scalarMap.to_rgba(color)
            if shaper == "C":
                circle=plt.Circle((x[ii]+size/2.0,y[ii]+size/2.0), size,color=color_val,ec="black",lw=1)
                plt.gca().add_artist(circle)
            if shaper == "R":
                recta=plt.Rectangle((x[ii],y[ii]), size, size, color=color_val,ec="black",lw=1)
                plt.gca().add_artist(recta)
        for ii in range(0, nspaxels):
            jj=spaxels[ii]
            poin=plt.Circle((x[jj]+size/2.0,y[jj]+size/2.0), 0.25,color='blue')
            plt.gca().add_artist(poin)
        ax2 = plt.subplot(gs[1])
        cb1 = colorbar.ColorbarBase(ax2, cmap=jet, norm=cNorm, orientation='vertical')
        cb1.set_label('Flux [%]')
        spec=t
        if i == 0:
            med=np.median(spec)
            y_max_s=3*med
            y_min_s=-0.5*med
        else:
            med=np.median(spec)
            y_max_s_now=3*med
            y_min_s_now=-0.5*med
            if y_max_s_now > y_max_s:
                y_max_s=y_max_s_now
            if y_min_s_now < y_min_s:
                y_min_s=y_min_s_now
        ax1 = plt.subplot(gs[2])
        plt.axis([wmin,wmax,y_min_s,y_max_s])
        plt.xlabel("Wavelength",fontsize=14)
        plt.ylabel("Flux",fontsize=14)
        plt.title("",fontsize=15)                      
        for ii in range(0, i+1):
            t_now=pdl_out[:,ii]
            spec=t_now
            ax1.plot(wave,spec)
        plt.tight_layout()
        if dev == "null":
            plt.show()
            plt.close()
        else:
            plt.savefig(dev)
            plt.close()
hd["CRPIX1"]=crpix 
hd["CRVAL1"]=crval 
hd["CDELT1"]=cdelt 
if pt.exists(output) == False:
    wfits(output,pdl_out,hd)
else:
    sycall("rm "+output)
    wfits(output,pdl_out,hd)

hd_e["CRPIX1"]=crpix 
hd_e["CRVAL1"]=crval 
hd_e["CDELT1"]=cdelt 
if pt.exists(e_output) == False:
    wfits(e_output,pdl_out_e,hd_e)
else:
    sycall("rm "+e_output)
    wfits(e_output,pdl_out_e,hd_e)

hd_m["CRPIX1"]=crpix 
hd_m["CRVAL1"]=crval 
hd_m["CDELT1"]=cdelt 
if pt.exists(weight_output) == False:
    wfits(weight_output,pdl_out_mask,hd_m)
else:
    sycall("rm "+weight_output)
    wfits(weight_output,pdl_out_mask,hd_m)