#!/usr/bin/python
import sys
import numpy as np
from numpy import float_
from numpy import absolute as abs
from numpy import random as ran
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
matplotlib.use('Agg')

def sycall(comand):
    from subprocess import call
    line=comand.split(" ")
    fcomand=[]
    fcomand.extend(line)
    call(fcomand)

def wfits(name, data, hdr):
    if ptt.exists(name) == False:
        wfit(name,data,hdr)
    else:
        sycall("rm "+name)
        wfit(name,data,hdr)
        
sys.argv=filter(None,sys.argv)
if len(sys.argv) < 6:
    print "USE: flux_elines_cube.py CUBE.fits,e_CUBE.fits,N_MC,CONT_CUBE.fits emission_lines.txt output.fits VEL_map.fits FWHM_map.fits[sigma_inst]"
    sys.exit(0)
speed_of_light=299792.458
junk_input_cube=sys.argv[1]
n_mc=1
if "," in junk_input_cube:
    data=junk_input_cube.split(",")
    input_cube=data[0]
    input_e_cube=data[1]
    n_mc=int(data[2])
    input_cont_cube=data[3]
else:
    input_cube=junk_input_cube
    input_e_cube="none"
elines_list=sys.argv[2]
output=sys.argv[3]
guided_map=sys.argv[4]
guided_sigma=sys.argv[5]

# Reading cube
[pdl_input_cube, h]=gdata(input_cube, 0, header=True)
[nz,nx,ny]=pdl_input_cube.shape
crpix=h["CRPIX3"]
crval=h["CRVAL3"]
cdelt=h["CDELT3"]
pdl_input_cube[np.isnan(pdl_input_cube)]=0
if input_e_cube != "none":
    if "[1]" in input_e_cube:
        input_e_cube=input_e_cube.replace('[1]','')
        pdl_input_e_cube=gdata(input_e_cube, 1)
    else:
        pdl_input_e_cube=gdata(input_e_cube)
    statst=my.stats(pdl_input_e_cube)
    pdl_input_e_cube[np.where(pdl_input_e_cube <= -5*statst[2])]=np.nan
    pdl_input_e_cube[np.where(pdl_input_e_cube >= +5*statst[2])]=np.nan
    pdl_input_e_cube[np.isnan(pdl_input_e_cube)]=5*statst[2]

if input_cont_cube != "none":
    pdl_input_cont_cube=gdata(input_cont_cube)
    statst=my.stats(pdl_input_cont_cube)
    pdl_input_cont_cube[np.isnan(pdl_input_cont_cube)]=statst[0]
#
# Reading velocity map
#
gmap=gdata(guided_map)
[nx1,ny1]=gmap.shape
if "fit" in guided_sigma:
    s_gmap=gdata(guided_sigma)
    s_gmap=s_gmap/2.354
else:
    s_gmap=float_(guided_sigma)*np.ones([nx1,ny1])
if nx1 < nx:
    nx=nx1
if ny1 < ny:
    ny=ny1
w_eline=[]
name=[]
ne=0
f=open(elines_list, 'r')
for line in f:
    if not "#" in line:
        data=line.split(" ")
        data=filter(None,data)
        w_eline.extend(float_([data[0]]))
        name.extend([data[1].replace('\n','')])
        ne=ne+1
f.close()
print str(ne)+" emission lines"
NZ_out=ne*4*2
pdl_out=np.zeros([NZ_out,nx,ny])
for k in range(0, ne):
    nn=k
    head="NAME"+str(nn)
    h[head]="flux "+name[k]
    nn=ne+k
    head="NAME"+str(nn)
    h[head]="vel "+name[k]
    nn=2*ne+k
    head="NAME"+str(nn)
    h[head]="disp "+name[k]
    nn=3*ne+k
    head="NAME"+str(nn)
    h[head]="EW "+name[k]
    #
    #  Errors
    #
    nn=4*ne+k
    head="NAME"+str(nn)
    h[head]="e_flux "+name[k]
    nn=5*ne+k
    head="NAME"+str(nn)
    h[head]="e_vel "+name[k]
    nn=6*ne+k
    head="NAME"+str(nn)
    h[head]="e_disp "+name[k]
    nn=7*ne+k
    head="NAME"+str(nn)
    h[head]="e_EW "+name[k]
    for i in range(0, nx):
        for j in range(0, ny):
        # We need to improved with an interpolated version!
            vel=gmap[i,j]
            sigma=s_gmap[i,j]
            sigma_ini=s_gmap[i,j]
            I0=0
            I1=0
            I2=0
            vel_I1=0
            s_I0=0
            s_I1=0
            s_I2=0
            s_vel_I1=0
            a_I0=np.zeros(n_mc)
            a_I1=np.zeros(n_mc)
            a_I2=np.zeros(n_mc)
            a_vel_I1=np.zeros(n_mc)
            f=(1+vel/speed_of_light)
            start_w=w_eline[k]*f-2.354*sigma
            end_w=w_eline[k]*f+2.354*sigma
            start_i=int((start_w-crval)/cdelt)
            end_i=int((end_w-crval)/cdelt)
            d_w=(end_w-start_w)/4.
            do_loop=1
            if start_i < 0:
                do_loop=0
            if end_i < 0:
                do_loop=0
            if start_i >= nz:
                do_loop=0
            if end_i >= nz:
                do_loop=0
            if sigma > 0 and do_loop == 1:
                for i_mc in range(0, n_mc):
                    I0=0
                    I1=0
                    I2=0
                    vel_I1=0
                    f=(1+vel/speed_of_light)
                    start_w=w_eline[k]*f-2.354*sigma
                    end_w=w_eline[k]*f+2.354*sigma
                    start_i=int((start_w-crval)/cdelt)
                    end_i=int((end_w-crval)/cdelt)
                    d_w=(end_w-start_w)/4.0
                    if start_i < 0:
                        start_i=0
                    if end_i < 0:
                        end_i=0
                    if start_i >= nz:
                        start_i=nz-1
                    if end_i >= nz:
                        end_i=nz-1
                    if i_mc == 0:
                        pdl_rnd=np.zeros(end_i-start_i+1)
                    else:
                        pdl_rnd=0.5-1.0*ran.rand(end_i-start_i+1)
                    n_I0=0
                    s_WE=0
                    sum_val=0
                    for iz in range(start_i, end_i):
                        val=pdl_input_cube[iz,i,j]        
                        if input_e_cube != "none":
                            e_val=pdl_input_e_cube[iz,i,j]
                        else:
                            e_val=0
                        val=val+(pdl_rnd[n_I0])*e_val
                        val=val*abs(sigma)*np.sqrt(2*3.1416)
                        w=crval+iz*cdelt
                        WE=np.exp(-0.5*((w-w_eline[k]*f)/(sigma))**2.0)
                        I0=I0+WE*val*np.exp(0.5*((w-w_eline[k]*f)/(sigma))**2.0)
                        I1=I1+val*w
                        s_WE=s_WE+WE
                        sum_val=sum_val+val
                        n_I0=n_I0+1
                    if n_I0 != 0:
                        s_WE=s_WE/n_I0
                        I0=I0/n_I0/s_WE
                        if sum_val != 0:
                            I1=I1/sum_val
                            if I1 != 0 and I0 != 0:
                                vel_I1=(I1/w_eline[k]-1.0)*speed_of_light
                            s_WE=0
                            n_I0
                            start_w=w_eline[k]*f-sigma
                            end_w=w_eline[k]*f+sigma
                            start_i=int((start_w-crval)/cdelt)
                            end_i=int((end_w-crval)/cdelt)
                            for iz in range(start_i, end_i):
                                val=pdl_input_cube[iz,i,j]            
                                val=val*abs(sigma)*np.sqrt(2*3.1416)
                                f1=(1.0+vel_I1/speed_of_light)
                                if I0 > 0 and val > 0 and I0 > val:
                                    S_now=(w-w_eline[k]*f1)/np.sqrt(2.0*(np.log(I0)-np.log(val)))
                                    if S_now > sigma:
                                        WE=1
                                        I2=I2+S_now*WE
                                        s_WE=s_WE+WE
                                        n_I0=n_I0+1
                            if s_WE > 0:
                                I2=I2/s_WE/2.345
                            else:
                                I2=sigma
                    a_I0[i_mc]=I0
                    a_I1[i_mc]=I1
                    a_I2[i_mc]=I2
                    a_vel_I1[i_mc]=vel_I1
            if n_mc > 1:
                I0=np.mean(a_I0)
                I1=np.mean(a_I1)
                I2=2.354*np.mean(a_I2)
                vel_I1=np.mean(a_vel_I1)
                s_I0=2.354*np.std(a_I0)
                s_I1=2.354*np.std(a_I1)
                s_I2=2.354*np.std(a_I2)
                s_vel_I1=2.354*np.std(a_vel_I1)
            else:
                I0=a_I0[0]
                I1=a_I1[0]
                I2=a_I2[0]
                vel_I1=a_vel_I1[0];
                s_I0=0
                s_I1=0
                s_I2=0
                s_vel_I1=0
            s_I0=abs(s_I0)
            s_I1=abs(s_I1)
            s_I2=abs(s_I2)
            s_vel_I1=abs(s_vel_I1)
            #
            # Continuum error!
            #
            start_i_0=int((start_w-60-crval)/cdelt)
            end_i_0=int((start_w-30-crval)/cdelt)
            start_i_1=int((end_w+30-crval)/cdelt)
            end_i_1=int((end_w+60-crval)/cdelt)
            if start_i_0 > 0 and end_i_0 > 0 and start_i_1 < nz and end_i_1 < nz:
                pdl_cont_0=pdl_input_cube[start_i_0:end_i_0,i,j]
                pdl_cont_1=pdl_input_cube[start_i_1:end_i_1,i,j]
                stats_cont_0=my.stats(pdl_cont_0)
                stats_cont_1=my.stats(pdl_cont_1)
                val_cont=0.5*(stats_cont_0[0]+stats_cont_1[1])
                e_val_cont=np.sqrt(((stats_cont_0[1])**2.0+(stats_cont_1[1])**2.0)/2.0)
                s_I0=np.sqrt(s_I0**2.0+(e_val_cont*2.354*I2)**2.0)
                s_I2=np.sqrt(s_I2**2.0+(s_I0/(2.354*I2))**2.0)
                s_vel_I1=np.sqrt(s_vel_I1**2.0+((s_I2/5500.0)*300000.0)**2.0)
            start_i_0=int((start_w-60.0-crval)/cdelt)
            end_i_0=int((start_w-30.0-crval)/cdelt)
            start_i_1=int((end_w+30.0-crval)/cdelt)
            end_i_1=int((end_w+60.0-crval)/cdelt)
            if start_i_0 > 0 and end_i_0 > 0 and start_i_1 < nz and end_i_1 < nz:
                pdl_cont_0=pdl_input_cont_cube[start_i_0:end_i_0,i,j]
                pdl_cont_1=pdl_input_cont_cube[start_i_1:end_i_1,i,j]
                stats_cont_0=my.stats(pdl_cont_0)
                stats_cont_1=my.stats(pdl_cont_1)
                val_cont=0.5*(stats_cont_0[0]+stats_cont_1[1])
                e_val_cont=np.sqrt(((stats_cont_0[1])**2.0+(stats_cont_1[1])**2.0)/2.0)
                if val_cont != 0:
                    EW=(-1.0)*I0/abs(val_cont)
                    e_EW=abs(s_I0)/abs(val_cont)+(I0*abs(e_val_cont))/(val_cont**2.0)
                    if abs(EW) < e_EW:
                        EW=0
                    if EW > 0:
                        EW=0
                else:
                    EW=np.nan
                    e_EW=np.nan
            else:
                EW=np.nan
                e_EW=np.nan
            pdl_out[0*ne+k,i,j]=I0
            pdl_out[1*ne+k,i,j]=vel_I1
            pdl_out[2*ne+k,i,j]=I2
            pdl_out[3*ne+k,i,j]=EW
            pdl_out[4*ne+k,i,j]=s_I0
            pdl_out[5*ne+k,i,j]=s_vel_I1
            pdl_out[6*ne+k,i,j]=s_I2
            pdl_out[7*ne+k,i,j]=e_EW
    print str(k)+"/"+str(ne)+","+str(w_eline[k])+","+name[k]+" DONE"
wfits(output, pdl_out, h)
sys.exit(0)
for k in ragne(0, ne):
    i_eline=int(w_eline[k])
    out_map="test_"+str(i_eline)+"_f.fits"
    pdl_map=pdl_out[k,:,:]
    wfits(out_map, pdl_map, h)
    kk=k+ne
    out_map="test_"+str(i_eline)+"_v.fits"
    pdl_map=pdl_out[kk,:,:]
    wfits(out_map, pdl_map, h)
    kk=k+2*ne
    out_map="test_"+str(i_eline)+"_d.fits"
    pdl_map=pdl_out[kk,:,:]
    wfits(out_map, pdl_map ,h)
