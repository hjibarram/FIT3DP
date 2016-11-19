#!/usr/bin/python

import sys
import numpy as np
from scipy import linalg as li
from numpy import float_
from numpy import absolute as abs
from numpy import random as ran
import matplotlib
from scipy.signal.signaltools import convolve2d
from scipy.interpolate.interpolate import interp1d
import my as my

def linear_least_squares(a, b, residuals=False):
    """
    Return the least-squares solution to a linear matrix equation.
    Solves the equation `a x = b` by computing a vector `x` that
    minimizes the Euclidean 2-norm `|| b - a x ||^2`.  The equation may
    be under-, well-, or over- determined (i.e., the number of
    linearly independent rows of `a` can be less than, equal to, or
    greater than its number of linearly independent columns).  If `a`
    is square and of full rank, then `x` (but for round-off error) is
    the "exact" solution of the equation.
    Parameters
    ----------
    a : (M, N) array_like
        "Coefficient" matrix.
    b : (M,) array_like
        Ordinate or "dependent variable" values.
    residuals : bool
        Compute the residuals associated with the least-squares solution
    Returns
    -------
    x : (M,) ndarray
        Least-squares solution. The shape of `x` depends on the shape of
        `b`.
    residuals : int (Optional)
        Sums of residuals; squared Euclidean 2-norm for each column in
        ``b - a*x``.
    """
    #  Copyright (c) 2013 Alexandre Drouin. All rights reserved.
    #  From https://gist.github.com/aldro61/5889795
    from warnings import warn
#    from scipy.linalg.fblas import dgemm
    from scipy.linalg.blas import dgemm
#    if type(a) != np.ndarray or not a.flags['C_CONTIGUOUS']:
#        warn('Matrix a is not a C-contiguous numpy array. The solver will create a copy, which will result' + \
#             ' in increased memory usage.')
    a = np.asarray(a, order='c')
    i = dgemm(alpha=1.0, a=a.T, b=a.T, trans_b=True)
    x = np.linalg.solve(i, dgemm(alpha=1.0, a=a.T, b=b)).flatten()
    if residuals:
        return x, np.linalg.norm(np.dot(a, x) - b)
    else:
        return x

def print_time():
    import time
    [year,mon,mday,hour,min,sec,wday,yday,isdst] = time.localtime()
    print "# TIME "+str(sec)+" "+str(min)+" "+str(hour)+" "+str(mday)+" "+str(mon)+" "+str(year)+" "+str(wday)+" "+str(yday)+" "+str(isdst)
    sec_now=hour*3600.+min*60.+sec
    return sec_now

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

def linfit1d(pdl_flux,pdl_model_n,weigh=[1]):
    #weigh=1.0/weigh#np.array(weigh)
    weigh=[1]
    [nx,ny]=pdl_model_n.shape
    pdl_flux=np.array([pdl_flux])
    if nx < ny:
        pdl_model_n=np.transpose(pdl_model_n)
        #pdl_flux=np.transpose(pdl_flux)
        nx=ny
    pdl_flux_m=pdl_flux/np.mean(pdl_flux)
    A=pdl_model_n
    B=pdl_flux_m
    if len(weigh) == nx:#Definicion del peso
        #print np.sum(pdl_flux.T)
        weigh=np.diag(weigh)
        A=np.dot(weigh,A)
        B=np.dot(weigh,np.transpose(B))
    else:
        B=np.transpose(B)
#    print B.shape
#    print A.shape
#    pdl_flux_n=np.zeros(nx)
#    pdl_flux_n[:]=pdl_flux[:,0]
#    pdl_flux=pdl_flux_n
#    coeffs_0=np.linalg.lstsq(A,pdl_flux,rcond=-1000)[0]
#    print coeffs_0
    #coeffs_0=li.lstsq(A,pdl_flux)[0]
    #coeffs_t=linear_least_squares(A,pdl_flux)
    #coeffs_0=np.zeros([len(coeffs_t),1])
    #coeffs_0[:,0]=coeffs_t
    #print coeffs_0
    coeffs_0= np.dot(np.linalg.inv(np.dot(A.T, A)),np.dot(A.T, B))*np.mean(pdl_flux)
    #print "____________"
    #print np.dot(np.linalg.inv(np.dot(A.T, A)),np.dot(A.T, A))
    #print coeffs_0
    #print "____________"
    #sys.exit()
    pdl_model_0=np.dot(A,coeffs_0)
    return [pdl_model_0,coeffs_0]

def plot_results(plot,pdl_wave,pdl_output,output_name,title):
    wave_now=pdl_wave
    if plot > 0 :
        if plot == 1:
            dev_plot="null"
            import matplotlib.pyplot as plt
            import matplotlib.colors as colors
            import matplotlib.cm as cmx
        else:
            matplotlib.use('Agg')
            import matplotlib.pyplot as plt
            import matplotlib.colors as colors
            import matplotlib.cm as cmx    
        if plot > 1:
            if plot == 2:
                dev_plot=output_name+".pdf"
            if plot == 3:
                dev_plot=output_name+".png"
            if plot > 3:
                dev_plot=output_name+".jpg"
        stats0=np.mean(pdl_output[:][0])
        stats1=np.sqrt(np.sum((pdl_output[:][0]-stats0)**2)/(len(pdl_output[:][0])-1))
        y_min=-0.5*stats0-0.75*stats1
        y_max=2*stats0+6*stats1
        nx=len(wave_now)
        ny=len(pdl_output)
        min_wave=wave_now[0]#5000.0#
        max_wave=wave_now[nx-2]#6200.0#
        fig = plt.figure()
        plt.axis([min_wave,max_wave,y_min,y_max])
        plt.xlabel("Wavelength",fontsize=14)
        plt.ylabel("Flux",fontsize=14)
        plt.title(title,fontsize=10)
        jet = plt.get_cmap('jet') 
        cNorm  = colors.Normalize(vmin=0, vmax=ny-1)
        scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)
        for j in range(0, ny):
            flux=pdl_output[:][j]
            color_val = scalarMap.to_rgba(j)
            plt.plot(wave_now,flux,color=color_val,lw=0.5)
        if dev_plot == "null":
            plt.show()
            plt.close()
        else:
            plt.savefig(dev_plot)
            plt.close()

def plot_results_min_max(plot,pdl_wave,pdl_output,output_name,title,y_min,y_max):
    wave_now=pdl_wave
    if plot > 0 :
        if plot == 1:
            dev_plot="null"
            import matplotlib.pyplot as plt
            import matplotlib.colors as colors
            import matplotlib.cm as cmx
        else:
            matplotlib.use('Agg')
            import matplotlib.pyplot as plt
            import matplotlib.colors as colors
            import matplotlib.cm as cmx    
        if plot > 1:
            if plot == 2:
                dev_plot=output_name+".pdf"
            if plot == 3:
                dev_plot=output_name+".png"
            if plot > 3:
                dev_plot=output_name+".jpg"
        nx=len(wave_now)
        ny=len(pdl_output)
        min_wave=wave_now[0]
        max_wave=wave_now[nx-2]
        fig = plt.figure()
        plt.axis([min_wave,max_wave,y_min,y_max])
        plt.xlabel("Wavelength",fontsize=14)
        plt.ylabel("Flux",fontsize=14)
        plt.title(title,fontsize=15)
        jet = plt.get_cmap('jet') 
        cNorm  = colors.Normalize(vmin=0, vmax=ny-1)
        scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)
        for j in range(0, ny):
            flux=pdl_output[:][j]
            color_val = scalarMap.to_rgba(j)
            plt.plot(wave_now,flux,color=color_val)
        if dev_plot == "null":
            plt.show()
            plt.close()
        else:
            plt.savefig(dev_plot)  
            plt.close()          

def smooth_ratio(pdl_data,pdl_model_spec_min,sigma):
    pdl_model_spec_min[np.where(pdl_model_spec_min == 0)[0]]=1.
    pdl_rat=pdl_data/pdl_model_spec_min
    pdl_rat[np.isnan(pdl_rat)]=1
    if sigma < 1 :
        sigma=1
    smooth_ratio=my.median_filter(int(7*2.354*sigma),pdl_rat)
    smooth_ratio[np.isnan(smooth_ratio)]=1
    return smooth_ratio

def fit_ssp_lin_no_zero(redshift,sigma,Av_NOW,crval,cdelt,crpix,nf,n_c,pdl_flux_c_ini,hdr_c_ini,wave_unc,masked,e_flux_unc,flux_unc,n_mc,chi_sq_min_now,MIN_CHISQ,plot):
    Av=Av_NOW
    iter_max=5
    last_chi=1e12
    n_unc=len(flux_unc)
    pdl_model_spec_min=np.zeros(n_unc)
    flux_masked=np.zeros(n_unc)
    coeffs=np.zeros([nf,3])    
    coeffs_N=np.zeros([nf,1])    
    coeffs_NM=np.zeros([nf,1])    
    coeffs_cat=np.zeros([nf+1,n_mc*n_mc])
    pdl_model_spec_cat=np.zeros([n_unc,n_mc*n_mc])

    for i in range(0, n_unc):
        flux_masked[i]=flux_unc[i]*masked[i]

    for i in range(0, nf):
        if Av[i] < 0:
            Av[i]=0

    wave_c=[]
    dpix_c_val=[]
    for j in range(0, n_c):
        wave_c.extend([(crval+cdelt*(j+1-crpix))*(1+redshift)])
        if j > 0:
            dpix_c_val.extend([wave_c[j]-wave_c[j-1]])


    dpix_c_val[0]=dpix_c_val[1]
    dpix_c=wave_c[1]-wave_c[0]
    rsigma=sigma/dpix_c

    name=[]
    age_mod=[]
    met_mod=[]
    ml=[]
    flux_c=[]
    model=np.zeros([n_unc,nf])
    model_no_mask=np.zeros([n_unc,nf])    
    age_min=0
    met_min=0
    med_norm=1
    Av_min=0
    age_min_mass=0
    met_min_mass=0
    Av_min_mass=0

    for iii in range(0, nf):
        header="NAME"+str(iii)
        name.extend([hdr_c_ini[header]]);
        name_min=name[iii]
        name_min=name_min.replace('spec_ssp_','')
        name_min=name_min.replace('.spec','')    
        name_min=name_min.replace('.dat','')
        data=name_min.split('_')
        AGE=data[0]
        MET=data[1]
        if 'Myr' in AGE:
            age=AGE.replace('Myr','')
            age=float_(age)/1000.
        else:
            age=AGE.replace('Gyr','')
            age=float_(age)
        met=float_(MET.replace('z','0.'))    
        age_mod.extend([age])
        met_mod.extend([met])
        header="NORM"+str(iii)    
        val_ml=float_(hdr_c_ini[header])
        if val_ml != 0:
            ml.extend([1/val_ml])
        else:
            ml.extend([1])
        box=int(3*rsigma)
        if box < 3:
            box=3
        kernel=np.zeros([1,2*box+1])
        norm=0
        flux_c.extend([0])
        for j in range(0, 2*box+1):
            gaus=np.exp(-0.5*(((j-box)/rsigma)**2))    
            kernel[0,j]=gaus
            norm=norm+gaus
        kernel=kernel/norm;
        pdl_flux_c_conv = convolve2d(pdl_flux_c_ini,kernel,mode='same')
        pdl_flux_c = pdl_flux_c_conv[iii,:]
        out_spec_pdl = interp1d(wave_c, pdl_flux_c,bounds_error=False,fill_value=0.)(wave_unc)
        n_c_out=out_spec_pdl.shape
        error=[]
        for i in range(0,n_unc):
            val=out_spec_pdl[i]
            if np.isnan(val):
                val=0
            model[i][iii]=val*masked[i]
            model_no_mask[i][iii]=val
            if masked[i] > 0:
                error.extend([0.01*abs(e_flux_unc[i])])
            else:
                error.extend([0])
    pdl_model=np.zeros([n_unc,nf])
    pdl_model_no_mask=np.zeros([n_unc,nf])
    pdl_error=np.zeros(n_unc)
    pdl_masked=masked
    pdl_dust_spec=np.zeros([n_unc,nf])
    for j in range(0, nf):
        for i in range(0, n_unc):
            wave_res=wave_unc[i]/(1+redshift)
            dust_rat=my.A_l(3.1,wave_res)
            dust=10**(-0.4*Av[j]*dust_rat)  
            pdl_dust_spec[i][j]=dust
            val=model[i][j]*dust
            pdl_model[i][j]=val
            val_no_mask=model_no_mask[i][j]*dust
            pdl_model_no_mask[i][j]=val_no_mask
            e_val=error[i]
            val_now=e_flux_unc[i]
            if val_now == 0:
                val_now=1
            pdl_error[i]=1.0/(abs(val_now)**2)
#
# We fit
#
    pdl_flux_masked=flux_masked
    ini_cat=0
    pdl_C_flast=np.ones(nf)

# Just a linear FIT, without restrictions!
    [y_model_now, coeffs] = linfit1d(pdl_flux_masked,pdl_model,1.0/pdl_error)
    y_model_now=y_model_now[:,0]
#
# We remove the models that are negative
#
    nf_new=0
    nf_neg=0
    for k in range(0, nf):
        C=coeffs[k][0]
        if C > 0:
            nf_new=nf_new+1
        else:
            nf_neg=nf_neg+1
    MOD=[]
    if nf_new > 0 :
        while nf_neg > 0:
            pdl_model_new=np.zeros([n_unc,nf_new])
            nf_i=0;
            for k in range(0, nf):
                C=coeffs[k][0]
                if C > 0:
                    pdl_model_new[:,nf_i] = pdl_model[:,k]
                    MOD.extend([k])
                    nf_i=nf_i+1
                else:
                    coeffs[k][0]=0
            [yfit, coeffs_new] = linfit1d(pdl_flux_masked,pdl_model_new,1.0/pdl_error)
            y_model_now=yfit[:,0]
            nf_i=0
            nf_neg=0
            nf_new=0
            for k in range(0, nf):
                C=coeffs[k][0]
                if C > 0:
                    val=coeffs_new[nf_i][0]
                    nf_i=nf_i+1
                    if val > 0:
                        coeffs[k][0]=val
                        nf_new=nf_new+1
                    else:
                        coeffs[k][0]=0
                        nf_neg=nf_neg+1
            if nf_new == 0:
                nf_neg=0
        #for k in range(0, nf):
        #    C=coeffs[k][0]
    else:
        nf_new=nf
##############################
# CHISQ VALUE
    chi=0
    chi2=0
    NFREE=0
    out_spec=[]
    chi_sec=[]
    res_spec=[]
    model_spec_min=[]
    model_spec=[]
    for j in range(0, n_unc):
        out_spec.extend([y_model_now[j]])
        model_spec.extend([out_spec[j]])
        res_spec.extend([flux_unc[j]-model_spec[j]])        
        model_spec_min.extend([model_spec[j]])
        chi_sec.extend([0])
        if flux_unc[j] != 0  and out_spec[j] !=0 and e_flux_unc[j] != 0:
            chi=chi+masked[j]*((flux_masked[j]-out_spec[j])**2.)/(e_flux_unc[j])**2.
            chi_sec[j]=masked[j]*((flux_unc[j]-out_spec[j])**2.)/(e_flux_unc[j])**2.
            NFREE=NFREE+1
    chi_sq=chi
    if NFREE > 0:
        chi_sq=(chi_sq/float_(NFREE))**0.5
    chi_sq_min_now=chi_sq
    y_model_end=y_model_now
    min_coeffs=coeffs
    chi_sq=chi_sq_min_now    
    norm_C=0
    norm_C_mass=0
    for k in range(0, nf):
        dust=10**(-0.4*Av[k]*dust_rat)
        C=coeffs[k][0]
        norm_C=norm_C+C
        norm_C_mass=norm_C_mass+C*ml[k]
        coeffs_N[k][0]=norm_C
        coeffs_NM[k][0]=norm_C_mass
        pdl_model_spec_min=pdl_model_spec_min+C*pdl_model_no_mask[:,k]    
    for k in range(0, nf):
        C=coeffs[k][0]
        if norm_C > 0:
            age_min=age_min+C*np.log10(age_mod[k])/norm_C
            met_min=met_min+C*met_mod[k]/norm_C
            Av_min=Av_min+C*Av[k]/norm_C        
            CN=C/norm_C
        C_now=C*med_norm
        if norm_C_mass > 0:
            age_min_mass=age_min_mass+C*np.log10(ml[k]*age_mod[k])/norm_C_mass
            met_min_mass=met_min_mass+C*ml[k]*met_mod[k]/norm_C_mass
            Av_min_mass=Av_min_mass+C*ml[k]*Av[k]/norm_C_mass
    age_min=10**(age_min)
    age_min_mass=10**(age_min_mass)
    
    pdl_age_mod=np.array(age_mod)
    pdl_met_mod=np.array(met_mod)
    pdl_ml=np.array(ml)
    pdl_Av=np.array(Av)

    pdl_res=flux_unc-pdl_model_spec_min
    pdl_wave_unc=wave_unc
    out_ps_now="junk";
    title="X="+str(chi_sq)+" Av="+str(Av[0])+" z="+str(redshift)+" sigma="+str(sigma)
    if plot > 0:
        plot_results(plot,pdl_wave_unc,[pdl_flux_masked,pdl_model_spec_min,pdl_res],out_ps_now,title)
    return [chi_sq,pdl_age_mod,pdl_met_mod,pdl_ml,pdl_Av,coeffs,coeffs_N,coeffs_NM,pdl_model_spec_min,pdl_res]

def fit_ssp_lin_no_zero_no_cont(redshift,sigma,Av_NOW,crval,cdelt,crpix,nf,n_c,pdl_flux_c_ini,hdr_c_ini,wave_unc,masked,e_flux_unc,flux_unc,n_mc,chi_sq_min_now,MIN_CHISQ,plot):
    Av=Av_NOW
    iter_max=5;
    last_chi=1e12
    n_unc=len(flux_unc)
    pdl_model_spec_min=np.zeros(n_unc)
    coeffs=np.zeros([nf,3])    
    coeffs_N=np.zeros([nf,1])    
    coeffs_NM=np.zeros([nf,1])    
    coeffs_cat=np.zeros([nf+1,n_mc*n_mc])
    pdl_model_spec_cat=np.zeros([n_unc,n_mc*n_mc])
    flux_masked=np.zeros(n_unc)
    for i in range(0, n_unc):
        flux_masked[i]=flux_unc[i]*masked[i]
    for i in range(0, nf):
        if Av[i] < 0:
            Av[i]=0
    wave_c=np.zeros(n_c)
    dpix_c_val=np.zeros(n_c)
    for j in range(0, n_c):
        wave_c[j]=(crval+cdelt*(j+1-crpix))*(1+redshift)
        if j > 0:
            dpix_c_val[j]=wave_c[j]-wave_c[j-1]
    dpix_c_val[0]=dpix_c_val[1]
    dpix_c=wave_c[1]-wave_c[0]
    rsigma=sigma/dpix_c
    name=[]
    age_mod=[]
    met_mod=[]
    ml=[]
    flux_c=[]
    model=np.zeros([n_unc,nf])
    model_no_mask=np.zeros([n_unc,nf])    
    age_min=0
    met_min=0
    Av_min=0
    age_min_mass=0
    met_min_mass=0
    Av_min_mass=0
    med_norm=1
    for iii in range(0, nf):
        header="NAME"+str(iii)
        name.extend([hdr_c_ini[header]]);
        name_min=name[iii]
        name_min=name_min.replace('spec_ssp_','')
        name_min=name_min.replace('.spec','')    
        name_min=name_min.replace('.dat','')
        data=name_min.split('_')
        AGE=data[0]
        MET=data[1]
        if 'Myr' in AGE:
            age=AGE.replace('Myr','')
            age=float_(age)/1000.
        else:
            age=AGE.replace('Gyr','')
            age=float_(age)
        met=float_(MET.replace('z','0.'))    
        age_mod.extend([age])
        met_mod.extend([met])
        header="NORM"+str(iii)    
        val_ml=float_(hdr_c_ini[header])
        if val_ml != 0:
            ml.extend([1/val_ml])
        else:
            ml.extend([1])
        box=int(3*rsigma)
        if box < 3:
            box=3
        kernel=np.zeros([1,2*box+1])
        norm=0    
        flux_c.extend([0])
        for j in range(0, 2*box+1):
            gaus=np.exp(-0.5*(((j-box)/rsigma)**2))    
            kernel[0,j]=gaus
            norm=norm+gaus
        kernel=kernel/norm;
        pdl_flux_c_conv = convolve2d(pdl_flux_c_ini,kernel,mode='same')
        pdl_flux_c = pdl_flux_c_conv[iii,:]
        out_spec_pdl = interp1d(wave_c, pdl_flux_c,bounds_error=False,fill_value=0.)(wave_unc)
        n_c_out=out_spec_pdl.shape
        error=[]
        for i in range(0,n_unc):
            val=out_spec_pdl[i]
            if np.isnan(val):
                val=0
            model[i][iii]=val*masked[i]
            model_no_mask[i][iii]=val
            if masked[i] > 0:
                error.extend([0.01*abs(e_flux_unc[i])])
            else:
                error.extend([0])
    pdl_model=np.zeros([n_unc,nf])
    pdl_model_no_mask=np.zeros([n_unc,nf])
    pdl_error=np.zeros(n_unc)
    pdl_masked=masked
    pdl_dust_spec=np.zeros([n_unc,nf])
    for j in range(0, nf):
        for i in range(0, n_unc):
            wave_res=wave_unc[i]/(1+redshift)
            dust_rat=my.A_l(3.1,wave_res)
            dust=10**(-0.4*Av[j]*dust_rat)  
            pdl_dust_spec[i][j]=dust
            val=model[i][j]*dust
            pdl_model[i][j]=val
            val_no_mask=model_no_mask[i][j]*dust
            pdl_model_no_mask[i][j]=val_no_mask
            e_val=error[i]
            val_now=e_flux_unc[i]
            if val_now == 0:
                val_now=1
            pdl_error[i]=1.0/(abs(val_now)**2)
#
# We fit
#
    pdl_flux_masked=flux_masked
    ini_cat=0
    pdl_C_flast=np.ones(nf)
# Just a linear FIT, without restrictions!
    [y_model_now, coeffs] = linfit1d(pdl_flux_masked,pdl_model,1/pdl_error)
    smooth_rati=smooth_ratio(flux_unc,y_model_now[:,0],int(sigma))
    y_model_now=y_model_now[:,0]*smooth_rati
#
# We remove the models that are negative
#
    nf_new=0
    nf_neg=0
    for k in range(0, nf):
        C=coeffs[k][0]
        if C > 0:
            nf_new=nf_new+1
        else:
            nf_neg=nf_neg+1
    MOD=[]
    if nf_new > 0 :
        while nf_neg > 0:
            pdl_model_new=np.zeros([n_unc,nf_new])
            nf_i=0;
            for k in range(0, nf):
                C=coeffs[k][0]
                if C > 0:
                    pdl_model_new[:,nf_i] = pdl_model[:,k]
                    MOD.extend([k])
                    nf_i=nf_i+1
                else:
                    coeffs[k][0]=0
            [yfit, coeffs_new] = linfit1d(pdl_flux_masked,pdl_model_new,1.0/pdl_error)
            smooth_rati=smooth_ratio(flux_unc,yfit[:,0],int(sigma))
            yfit=yfit[:,0]*smooth_rati
            y_model_now=yfit
            nf_i=0
            nf_neg=0
            nf_new=0
            for k in range(0, nf):
                C=coeffs[k][0]
                if C > 0:
                    val=coeffs_new[nf_i][0]
                    nf_i=nf_i+1
                    if val > 0:
                        coeffs[k][0]=val
                        nf_new=nf_new+1
                    else:
                        coeffs[k][0]=0
                        nf_neg=nf_neg+1
            if nf_new == 0:
                nf_neg=0
        #for k in range(0, nf):
        #    C=coeffs[k][0]
    else:
        nf_new=nf
    #sys.exit()
##############################
# CHISQ VALUE
    chi=0
    chi2=0
    NFREE=0
    out_spec=[]
    chi_sec=[]
    res_spec=[]
    model_spec_min=[]
    model_spec=[]
    for j in range(0, n_unc):
        out_spec.extend([y_model_now[j]])
        model_spec.extend([out_spec[j]])
        res_spec.extend([flux_unc[j]-model_spec[j]])        
        model_spec_min.extend([model_spec[j]])
        chi_sec.extend([0])
        if flux_unc[j] != 0  and out_spec[j] !=0 and e_flux_unc[j] != 0:
            chi=chi+masked[j]*((flux_masked[j]-out_spec[j])**2.)/(e_flux_unc[j])**2.
            chi_sec[j]=masked[j]*((flux_unc[j]-out_spec[j])**2.)/(e_flux_unc[j])**2.
            NFREE=NFREE+1
    chi_sq=chi
    if NFREE > 0:
        chi_sq=(chi_sq/float_(NFREE))**0.5
    chi_sq_min_now=chi_sq
    y_model_end=y_model_now
    min_coeffs=coeffs
    chi_sq=chi_sq_min_now    
    norm_C=0
    norm_C_mass=0
    for k in range(0, nf):
        dust=10**(-0.4*Av[k]*dust_rat)
        C=coeffs[k][0]
        norm_C=norm_C+C
        norm_C_mass=norm_C_mass+C*ml[k]
        coeffs_N[k][0]=norm_C
        coeffs_NM[k][0]=norm_C_mass
        pdl_model_spec_min=pdl_model_spec_min+C*pdl_model_no_mask[:,k]  
    smooth_rati=smooth_ratio(flux_unc,pdl_model_spec_min,int(sigma))
    pdl_model_spec_min=pdl_model_spec_min*smooth_rati
    for k in range(0, nf):
        C=coeffs[k][0]
        if norm_C > 0:
            age_min=age_min+C*np.log10(age_mod[k])/norm_C
            met_min=met_min+C*met_mod[k]/norm_C
            Av_min=Av_min+C*Av[k]/norm_C        
            CN=C/norm_C
        C_now=C*med_norm
        if norm_C_mass > 0:
            age_min_mass=age_min_mass+C*np.log10(ml[k]*age_mod[k])/norm_C_mass
            met_min_mass=met_min_mass+C*ml[k]*met_mod[k]/norm_C_mass
            Av_min_mass=Av_min_mass+C*ml[k]*Av[k]/norm_C_mass
    age_min=10**(age_min)
    age_min_mass=10**(age_min_mass)
    
    pdl_age_mod=np.array(age_mod)
    pdl_met_mod=np.array(met_mod)
    pdl_ml=np.array(ml)
    pdl_Av=np.array(Av)

    pdl_res=flux_unc-pdl_model_spec_min
    pdl_wave_unc=wave_unc
    out_ps_now="junk";
    title="X="+str(chi_sq)+" Av="+str(Av[0])+" z="+str(redshift)+" sigma="+str(sigma)
    if plot > 0:
#        plot_results(plot,pdl_wave_unc,[pdl_flux_masked,pdl_model_spec_min,pdl_res],out_ps_now,title)
        plot_results(plot,pdl_wave_unc,[pdl_flux_masked,y_model_end,pdl_res],out_ps_now,title)
        #sys.exit()
    return [chi_sq,pdl_age_mod,pdl_met_mod,pdl_ml,pdl_Av,coeffs,coeffs_N,coeffs_NM,pdl_model_spec_min,pdl_res]

def copy_a(n_mod,a_print):
    a_copy=np.zeros([n_mod, 9])
    for ii in range(0, n_mod):   
        for j in range(0, 9):
            a_copy[ii][j]=a_print[ii][j];
    return a_copy

def copy_a_results(n_mod,kk,a_print,a_copy):
    if kk > 0:
        a_copy=np.concatenate((a_copy,np.zeros([1, n_mod, 9])),axis=0)
    for ii in range(0, n_mod):   
        for j in range(0, 9):
            a_copy[kk][ii][j]=a_print[ii][j]
    return a_copy    

def create_single_model(pdl_wave,i_now,typef,a_c):
    dims=pdl_wave.shape
    nx=dims[0]
    pdl_out=np.zeros([1,nx])
    for i in range(0, nx):
        w=pdl_wave[i]
        if typef[i_now] == "eline\n":
            speed_of_light=299792.458
            factor=(1+a_c[i_now][3]/speed_of_light+a_c[i_now][5])
            e1=1
            Y1=1
            if a_c[i_now][2] != 0:
                e1=np.exp(-0.5*((w-a_c[i_now][0]*factor)/a_c[i_now][2])**2)
                Y1=a_c[i_now][1]*e1/(a_c[i_now][2]*((2*3.1416)**0.5))
            pdl_out[0][i]=Y1
        if typef[i_now] == "poly1d\n":
            Yi=0
            for ii in range(0, 9):
                Yi=Yi+a_c[i_now][ii]*(w)**(ii)
            pdl_out[0][i]=Yi
    return pdl_out

def create_single_model_one(pdl_wave,i_now,typef,a_c):
    dims=pdl_wave.shape
    nx=dims[0]
    pdl_out=np.zeros([1,nx])
    for i in range(0, nx):
        w=pdl_wave[i]
        if typef[i_now] == "eline\n":
            speed_of_light=299792.458
            factor=(1+a_c[i_now][3]/speed_of_light+a_c[i_now][5])
            e1=1
            Y1=1
            if a_c[i_now][2] != 0:
                e1=np.exp(-0.5*((w-a_c[i_now][0]*factor)/a_c[i_now][2])**2)
                Y1=1.0*e1/(a_c[i_now][2]*((2*3.1416)**0.5))
            pdl_out[0][i]=Y1
        if typef[i_now] == "poly1d\n":
            Yi=0
            for ii in range(0, 9):
                Yi=Yi+a_c[i_now][ii]*(w)**(ii)
            pdl_out[0][i]=Yi
    return pdl_out

def create_single_model_poly(pdl_wave,ii):
    dims=pdl_wave.shape
    nx=dims[0]
    pdl_out=np.zeros([1,nx])
    for i in range(0, nx):
        w=pdl_wave[i]
        Yi=(w)**(ii)
        pdl_out[0][i]=Yi
    return pdl_out
 
def fit_elines_grad_rnd_new(pdl_wave,pdl_flux,pdl_e_flux,n_mod,chi_goal,d_chi_goal,typef,a,ia,a0,a1,link,n_rnd,pdl_masked,defi,SCALE_INI):
    ttt=1
    a_out=np.zeros([n_mod, 9])    
    SCALE=SCALE_INI
    n_mc=n_rnd
    cont=0
    dims=pdl_flux.shape    
    if len(dims) == 1:
        nx=dims[0]

    # First Guess
    n_points=np.sum(pdl_masked)
    pdl_model_0=np.zeros(nx)
    pdl_model_cont_0=np.zeros(nx)
    for i in range(0,n_mod):
        pdl_tmp=create_single_model(pdl_wave,i,typef,a)
        pdl_model_0=pdl_model_0+pdl_tmp
    
    pdl_model_cont_0=create_single_model(pdl_wave,(n_mod-1),typef,a)
    pdl_chi_now=pdl_masked*((pdl_flux-pdl_model_0)**2)/(pdl_e_flux**2)
    chi_sq=np.sum(pdl_chi_now)
    chi_sq=3*(chi_sq/(n_points-n_mod-1))**0.5
    chi_sq_ini=chi_sq
    a_out=copy_a(n_mod,a)
    a_now=np.zeros([n_mod, 9])
    pdl_rnd=ran.randn(n_mod*9*n_mc)
    pdl_rnd_lin=ran.rand(n_mod*9*n_mc)
# 1st we derive the redshift!
    new_mc=int(n_mc/2)
    for ii in range(0, new_mc):
    # We change slighly the parameters
    #
        for i in range(0, n_mod):
            for j in range(0, 9):
                rnd=pdl_rnd[(j+(9*i)+(9*n_mod*ii))]#*ii]
                rnd_lin=pdl_rnd_lin[(j+(9*i)+(9*n_mod*ii))]#*ii]
                if ia[i][j]==1:
                    if link[i][j]==-1:
                        if typef[i] != "poly1d\n":
                            if j==3:
                                a_now[i][j]=a0[i][j]+rnd_lin*(a1[i][j]-a0[i][j])#*ii/new_mc                                
                            else:
                                a_now[i][j]=a_out[i][j]                
                        else:
                            a_now[i][j]=a_out[i][j]+SCALE*a_out[i][j]*rnd
                        
                        if a_now[i][j] < a0[i][j]:
                            a_now[i][j] = a0[i][j]
                        if a_now[i][j] > a1[i][j]:
                            a_now[i][j] = a1[i][j]
                    else:
                        k=link[i][j]-1
                        method=a1[i][j]
                        if method==0:
                            a_now[i][j]=a_now[k][j]+a0[i][j]
                        else:
                            a_now[i][j]=a_now[k][j]*a0[i][j]
                else:
                    a_now[i][j]=a_out[i][j]                    
        n_mod_free=0
        for i1 in range(0,n_mod):    
            if typef[i1] == "eline\n":
                if ia[i1][1] == 1:
                    if link[i1][1] == -1:
                        pdl_tmp=create_single_model_one(pdl_wave,i1,typef,a_now)
                        if n_mod_free == 0:
                            pdl_model_n=pdl_tmp
                        else:
                            pdl_model_n=np.vstack((pdl_model_n,pdl_tmp))           
                        n_mod_free=n_mod_free+1
                    else:
                        pdl_tmp=create_single_model_one(pdl_wave,i1,typef,a_now)
                        ik=link[i1][1]-1
                        t = pdl_model_n[ik,:]
                        pdl_model_n[ik,:]=pdl_tmp*a0[i1][1]+t
            if typef[i1] == "poly1d\n":        
                for j in range(0, 9):
                    if ia[i1][j] == 1:
                        pdl_tmp=create_single_model_poly(pdl_wave,j)
                        pdl_model_n=np.vstack((pdl_model_n,pdl_tmp))
                        n_mod_free=n_mod_free+1
        pdl_grandom=ran.randn(nx)
        pdl_flux_fit=pdl_flux-cont+pdl_grandom*pdl_e_flux
        dim_n=pdl_model_n.shape
        if dim_n[0] > 1:
            [pdl_model_0, coeffs_0] = linfit1d(pdl_flux_fit,pdl_model_n)
        else:
            pdl_model_n[np.where(pdl_model_n == 0)[0]]=1.
            pdl_rat=pdl_flux_fit/pdl_model_n
            pdl_rat[np.isnan(pdl_rat)]=0
            stats=np.median(pdl_rat)
            coeffs_0=np.zeros(1)+1
            pdl_model_0=stats*pdl_model_n
            coeffs_0=stats*coeffs_0
        
        pdl_model_0=pdl_model_0+cont
        n_mod_free=0
        for i1 in range(0, n_mod):
            if typef[i1] == "eline\n":
                if ia[i1][1] == 1:
                    if link[i1][1] == -1:
                        a_now[i1][1] = coeffs_0[n_mod_free]
                        if a_now[i1][1] < a0[i1][1]:
                            a_now[i1][1]=a0[i1][1]
                        if a_now[i1][1] > a1[i1][1]:
                            a_now[i1][1]=a1[i1][1]
                        n_mod_free=n_mod_free+1
                    else:
                        k=link[i1][1]-1
                        method=a1[i1][1]
                        if method == 0:
                            a_now[i1][1]=a_now[k][1]+a0[i1][1]
                        else:
                            a_now[i1][1]=a_now[k][1]*a0[i1][1]
            if typef[i1] == "poly1d\n":
                for j in range(0, 9):
                    if ia[i1][j] == 1:
                        val=coeffs_0[n_mod_free]
                        a_now[i1][j]=val
                        n_mod_free=n_mod_free+1
        
        pdl_model_0=np.zeros(nx)
        pdl_model_cont_0=np.zeros(nx)
        #pdl_model=pdl_model_0
        #pdl_model_cont=pdl_model_cont_0
        for i in range(0, n_mod):
            pdl_tmp=create_single_model(pdl_wave,i,typef,a_now)
            pdl_model_0=pdl_model_0+pdl_tmp
        pdl_model_cont_0=create_single_model(pdl_wave,(n_mod-1),typef,a_now)
        pdl_a_now=a_now
    
        pdl_chi_now=pdl_masked*((pdl_flux-pdl_model_0)**2)/(pdl_e_flux**2)
        chi_sq=np.sum(pdl_chi_now)
        chi_sq=(chi_sq/(n_points-n_mod-1))**0.5
        
        if chi_sq <= chi_sq_ini:
            a_out=copy_a(n_mod,a_now)
            chi_sq_ini=chi_sq
            pdl_model=pdl_model_0
            pdl_model_cont=pdl_model_cont_0
            ttt=0
    ii=0
    a_now=copy_a(n_mod,a_out)
    chi_sq_ini=1e12;
    #print chi_sq, chi_sq_ini
#
# 2nd we derive the sigma
#
    new_mc=int(n_mc/3)
#    pdl_rnd=ran.randn(n_mod*9*n_mc)
#    pdl_rnd_lin=ran.rand(n_mod*9*n_mc)

    for ii in range(0, new_mc):
        #
        # We change slighly the parameters
        #
        for i in range(0, n_mod):
            for j in range(0, 9):
                rnd=pdl_rnd[(j+(9*i)+(9*n_mod*ii))]
                rnd_lin=pdl_rnd_lin[(j+(9*i)+(9*n_mod*ii))]
                if ia[i][j] == 1:
                    if link[i][j] == -1:
                        if typef[i] != "poly1d\n":
                            if j == 2:
                                a_now[i][j]=a0[i][j]+rnd_lin*(a1[i][j]-a0[i][j])#*ii/new_mc
                                #print a_now[i][j], a0[i][j], (a1[i][j]-a0[i][j]), chi_sq
                            else:
                                a_now[i][j]=a_out[i][j]   
                        else:
                            a_now[i][j] = a_out[i][j]+SCALE*a_out[i][j]*rnd
                        if a_now[i][j] < a0[i][j]:
                            a_now[i][j] = a0[i][j]
                        if a_now[i][j] > a1[i][j]:
                            a_now[i][j] = a1[i][j]
                    else:
                        k=link[i][j]-1
                        method=a1[i][j]
                        if method == 0:
                            a_now[i][j] = a_now[k][j]+a0[i][j]
                        else:
                            a_now[i][j] = a_now[k][j]*a0[i][j]
                else:
                    a_now[i][j]=a_out[i][j]

        n_mod_free=0
        for i1 in range(0, n_mod):
            if typef[i1] == "eline\n":
                if ia[i1][1] == 1:
                    if link[i1][1] == -1:
                        pdl_tmp=create_single_model_one(pdl_wave,i1,typef,a_now)
                        if n_mod_free == 0:
                            pdl_model_n=pdl_tmp
                        else:
                            pdl_model_n=np.vstack((pdl_model_n,pdl_tmp))
                        n_mod_free=n_mod_free+1
                    else:
                        pdl_tmp=create_single_model_one(pdl_wave,i1,typef,a_now)
                        i2=link[i1][1]-1
                        t=pdl_model_n[i2,:]
                        pdl_model_n[i2,:]=pdl_tmp*a0[i1][1]+t
            if typef[i1] == "poly1d\n":
                for j in range(0, 9):
                    if ia[i1][j] == 1:
                        pdl_tmp=create_single_model_poly(pdl_wave,j)
                        pdl_model_n=np.vstack((pdl_model_n,pdl_tmp))
                        n_mod_free=n_mod_free+1
                        
        pdl_grandom=ran.randn(nx)
        pdl_flux_fit=pdl_flux-cont+pdl_grandom*pdl_e_flux
        dim_n=pdl_model_n.shape
        if dim_n[0] > 1:
            [pdl_model_0, coeffs_0] = linfit1d(pdl_flux_fit,pdl_model_n)
        else:
            pdl_model_n[np.where(pdl_model_n == 0)[0]]=1.
            pdl_rat=pdl_flux_fit/pdl_model_n
            pdl_rat[np.isnan(pdl_rat)]=0
            stats=np.median(pdl_rat)
            coeffs_0=np.zeros(1)+1
            pdl_model_0=stats*pdl_model_n
            coeffs_0=stats*coeffs_0
        pdl_model_0=pdl_model_0+cont
        n_mod_free=0
        for i1 in range(0, n_mod):
            if typef[i1] == "eline\n":
                if ia[i1][1] == 1:
                    if link[i1][1] == -1:
                        a_now[i1][1]=coeffs_0[n_mod_free]
                        n_mod_free=n_mod_free+1
                    else:
                        k=link[i1][1]-1
                        method=a1[i1][1]
                        if method == 0:
                            a_now[i1][1]=a_now[k][1]+a0[i1][1]
                        else:
                            a_now[i1][1]=a_now[k][1]*a0[i1][1]
                for j in range(0, 9):
                    if ia[i1][j] == 1 and link[i1][j] == -1:
                        if a_now[i1][j] < a0[i1][j]:
                            a_now[i1][j]=a0[i1][j]
                        if a_now[i1][j] > a1[i1][j]:
                            a_now[i1][j]=a1[i1][j]
            if typef[i1] == "poly1d\n":
                for j in range(0, 9):
                    if ia[i1][j] == 1:
                        val=coeffs_0[n_mod_free]
                        a_now[i1][j]=val
                        n_mod_free=n_mod_free+1    
        pdl_model_cont_0=create_single_model(pdl_wave,(n_mod-1),typef,a_now)
        pdl_a_now=a_now
        #print pdl_model_cont_0.shape
        pdl_chi_now=pdl_masked*((pdl_flux-pdl_model_0[:,0])**2)/(pdl_e_flux**2)
        chi_sq=np.sum(pdl_chi_now)
        chi_sq=(chi_sq/(n_points-n_mod-1))**0.5
        if chi_sq <= chi_sq_ini:
            a_out=copy_a(n_mod,a_now)
            chi_sq_ini=chi_sq
            pdl_model=pdl_model_0
            pdl_model_cont=pdl_model_cont_0
            ttt=0
            #print a_out[0][1],chi_sq
    a_now=copy_a(n_mod,a_out)
    # We start the fitting loop!

    ii=0
    SCALE_IN=SCALE
    new_mc=0
    while ii < new_mc:
    #
    # We change slighly the parameters
    #
        for i in range(0, n_mod):
            for j in range(0, 9):
                rnd=pdl_rnd[(j+(9*i))*ii]
                rnd_lin=pdl_rnd_lin[(j+(9*i))*ii]
                if ia[i][j] == 1:
                    if link[i][j] == -1:
                        A1=a1[i][j]
                        A0=a0[i][j]
                        if a1[i][j] > 1.3*a_out[i][j]:
                            A1=1.3*a_out[i][j]
                        if a0[i][j] < 0.7*a_out[i][j]:
                            A0=0.7*a_out[i][j]
                        if typef[i] == "eline\n":
                            if j == 3:
                                a_now[i][j]=a_out[i][j]+SCALE*rnd*(A1-A0)/(5*new_mc)
                            else:
                                a_now[i][j]=a_out[i][j]+SCALE*a_out[i][j]*rnd
                            if a_now[i][j] < A0:
                                a_now[i][j]=A0
                            if a_now[i][j] > A1:
                                a_now[i][j]=A1
                        else:
                            a_now[i][j]=a_out[i][j]
                    else:
                        k=link[i][j]-1
                        method=a1[i][j]
                        if method == 0:
                            a_now[i][j]=a_now[k][j]+a0[i][j]
                        else:
                            a_now[i][j]=a_now[k][j]*a0[i][j]
                else:
                    a_now[i][j]=a_out[i][j]
        # MC test
        pdl_model_0=np.zeros(nx)
        pdl_model_cont_0=np.zeros(nx)
        for i in range(0, n_mod):
            pdl_tmp=create_single_model(pdl_wave,i,typef,a_now)
            pdl_model_0=pdl_model_0+pdl_tmp
        pdl_model_cont_0=create_single_model(pdl_wave,(n_mod-1),typef,a_now)
        pdl_chi_now=pdl_masked*((pdl_flux-pdl_model_0)**2)/(pdl_e_flux**2)
        chi_sq=np.sum(pdl_chi_now)
        chi_sq=(chi_sq/(n_points-n_mod-1))**0.5
        if chi_sq < chi_sq_ini:
            a_out=copy_a(n_mod, a_now)
            chi_sq_ini=chi_sq
            pdl_model=pdl_model_0
            pdl_model_cont=pdl_model_cont_0
            SCALE=SCALE*0.99
            ttt=0
            if SCALE < 0.1*SCALE_IN:
                SCALE=SCALE_IN*0.1
        else:
            SCALE=SCALE_IN
            if (abs(chi_sq-chi_sq_ini) < d_chi_goal) or (chi_sq_ini < chi_goal): 
                ii=n_mc
        ii=ii+1   
    chi_sq_now=chi_sq_ini
    if ttt == 1:
        pdl_model=pdl_model_0
        pdl_model_cont=pdl_model_cont_0
    #print ttt, "HOLA  TEST!!!!!!!"
    #
    # We force lineal!
    #
    #a_out=copy_a(n_mod,a_now_lin)
    pdl_a_out=a_out
    return [chi_sq_now,pdl_a_out,pdl_model,pdl_model_cont]

def copy_a_pdl(n_mod,a_print):
    a_copy=np.zeros([n_mod, 9])
    for ii in range(0, n_mod):   
        for j in range(0, 9):
            a_copy[ii][j]=a_print[ii][j]
    return a_copy

def print_a_final(n_mod,a_print,typef,chi_a,back_noise=[1],si=0):
    nx=len(back_noise)
    if nx == 1:
        back_noise=0
    print str(n_mod)+" "+str(chi_a)
    line=""
    for ii in range(0, n_mod):
        if typef[ii] == "eline\n":
            typefs="eline"
        else:
            typefs="poly1d"
        line=line+typefs+" "
        e_F=back_noise*2.354*a_print[0][ii][2]
        a_print[1][ii][1]=np.sqrt((a_print[1][ii][1])**2+e_F**2)
        for j in range(0, 9):
            line=line+str(a_print[0][ii][j])+" "+str(a_print[1][ii][j])+" "
        if si == 0:
            print line
        line=""
    return a_print

def print_a_final_file(n_mod,a_print,typef,chi_a,outfile):
    print "OUT_FILE = "+outfile
    f=open(outfile, "w")
    f.write(str(n_mod)+" "+str(chi_a)+"\n")
    line=""
    for ii in range(0, n_mod):   
        if typef[ii] == "eline\n":
            typefs="eline"
        else:
            typefs="poly1d"
        line=line+typefs+" "
        for j in range(0, 9):
            line=line+str(a_print[0][ii][j])+" "+str(a_print[1][ii][j])+" "
        f.write(line+"\n")
        line=""
    f.close

def print_a_final_file_add(n_mod,a_print,typef,chi_a,outfile,i_val=0,sile=0):
    if sile == 0:
        print "OUT_FILE = "+outfile
    f=open(outfile, "a")
    f.write(str(n_mod)+" "+str(chi_a)+" "+str(i_val)+"\n")
    line=""
    for ii in range(0, n_mod):
        if typef[ii] == "eline\n":
            typefs="eline"
        else:
            typefs="poly1d"
        line=line+typefs+" "
        for j in range(0, 9):
            line=line+str(a_print[0][ii][j])+" "+str(a_print[1][ii][j])+" "
        f.write(line+"\n")
        line=""
    f.close

def print_a_final_file_add_mpi(n_mod,a_print,typef,chi_a,outfile,i_val=0):
    #print "OUT_FILE = "+outfile
    #f=open(outfile, "a")
    #f.write(str(n_mod)+" "+str(chi_a)+" "+str(i_val)+"\n")
    warr=[]
    warr.extend([str(n_mod)+" "+str(chi_a)+" "+str(i_val)+"\n"])
    line=""
    for ii in range(0, n_mod):
        if typef[ii] == "eline\n":
            typefs="eline"
        else:
            typefs="poly1d"
        line=line+typefs+" "
        for j in range(0, 9):
            line=line+str(a_print[0][ii][j])+" "+str(a_print[1][ii][j])+" "
        #f.write(line+"\n")
        warr.extend([line+"\n"])
        line=""
    #f.close
    return warr

def add_back_noise(n_mod,a_print,typef,chi_a,back_noise):
    for ii in range(0, n_mod):
        e_F=back_noise*2.354*a_print[0][ii][2]
        a_print[1][ii][1]=np.sqrt((a_print[1][ii][1])**2+e_F**2)
    return a_print

def mean_a_results(n_mod,nk,a_print,ia):
    a_copy=np.zeros([2,n_mod,9])
    for ii in range(0, n_mod):
        for j in range(0, 9):
            a_tmp=np.zeros(nk)
            for k in range(0, nk): 
                a_tmp[k]=a_print[k][ii][j]
            val=np.mean(a_tmp)
            e_val=np.std(a_tmp)
            a_copy[0][ii][j]=val
            a_copy[1][ii][j]=e_val*ia[ii][j]
    return a_copy

def mean_a_results_last(n_mod,nk,a_print,ia):
    if nk == 0:
        nk=1
    a_copy=np.zeros([2,n_mod,9])
    for ii in range(0, n_mod):   
        for j in range(0, 9):
            a_tmp=np.zeros(nk)
            for k in range(0, nk): 
                a_tmp[k]=a_print[k][ii][j]
            val=np.mean(a_tmp)
            e_val=np.std(a_tmp)
            if e_val == 0:
                e_val=0.1*val
            a_copy[0][ii][j]=a_tmp[nk-1]
            a_copy[1][ii][j]=e_val*ia[ii][j]
    return a_copy

def add_a_results_elines(n_mod,a_final,a_type,n_mod_fixed,a_final_fixed,a_type_fixed):
    KK=0
    for ii in range(0, n_mod):   
        if a_type[ii] == "eline\n":
            if n_mod_fixed+KK > 0:
                a_final_fixed=np.vstack((a_final_fixed,np.zeros(9)))
            a_type_fixed.extend([a_type[ii]])    
            for j in range(0, 9):
                a_final_fixed[KK+n_mod_fixed][j]=a_final[0][ii][j]
            KK=KK+1
    n_mod_fixed=n_mod_fixed+KK
    return [n_mod_fixed,a_final_fixed,a_type_fixed]

def fit_ssp_lin_MC(redshift,sigma,Av_NOW,crval,cdelt,crpix,nf,n_c,pdl_flux_c_ini,hdr_c_ini,wave_unc,masked,e_flux_unc,flux_unc,n_mc,chi_sq_min_now,MIN_CHISQ,plot):
    Av=Av_NOW;
    iter_max=5
    last_chi=1e12
    flux_unc=np.array(flux_unc)
    flux_unc[np.isnan(flux_unc)]=0
    flux_unc[np.isinf(flux_unc)]=0
    n_unc=len(flux_unc)
    coeffs_N=np.zeros([nf,1])    
    coeffs_NM=np.zeros([nf,1])    
    coeffs_cat=np.zeros([nf+1,n_mc])
    pdl_model_spec_cat=np.zeros([n_unc,n_mc])
    pdl_1st_model=[]
    flux_masked=np.zeros(n_unc)
    for i in range(0, n_unc):
        flux_masked[i]=flux_unc[i]*masked[i]
    for i in range(0, nf):
        if Av[i] < 0:
            Av[i]=0
    wave_c=np.zeros(n_c)
    dpix_c_val=[]
    for j in range(0, n_c):
        wave_c[j]=(crval+cdelt*(j+1-crpix))*(1+redshift)
        if j > 0:
            dpix_c_val.extend([wave_c[j]-wave_c[j-1]])    
    dpix_c_val[0]=dpix_c_val[1]
    dpix_c=wave_c[1]-wave_c[0]
    rsigma=sigma/dpix_c
    name=[]
    age_mod=[]
    met_mod=[]
    ml=[]
    flux_c=[]
    model=np.zeros([n_unc,nf])
    model_no_mask=np.zeros([n_unc,nf])   
    MED_NORM=[]
    age_min=0
    met_min=0
    Av_min=0
    age_min_mass=0
    met_min_mass=0
    Av_min_mass=0
    for iii in range(0, nf):
        header="NAME"+str(iii)
        name.extend([hdr_c_ini[header]]);
        name_min=name[iii]
        name_min=name_min.replace('spec_ssp_','')
        name_min=name_min.replace('.spec','')    
        name_min=name_min.replace('.dat','')
        data=name_min.split('_')
        AGE=data[0]
        MET=data[1]
        if 'Myr' in AGE:
            age=AGE.replace('Myr','')
            age=float_(age)/1000.
        else:
            age=AGE.replace('Gyr','')
            age=float_(age)
        met=float_(MET.replace('z','0.'))    
        age_mod.extend([age])
        met_mod.extend([met])
        header="NORM"+str(iii)    
        val_ml=float_(hdr_c_ini[header])
        if val_ml != 0:
            ml.extend([1/val_ml])
        else:
            ml.extend([1])
        box=int(3*rsigma)
        if box < 3:
            box=3
        kernel=np.zeros([1,2*box+1])
        norm=0    
        flux_c.extend([0])
        for j in range(0, 2*box+1):
            gaus=np.exp(-0.5*(((j-box)/rsigma)**2))    
            kernel[0,j]=gaus
            norm=norm+gaus
        kernel=kernel/norm;
        pdl_flux_c_conv = convolve2d(pdl_flux_c_ini,kernel,mode='same')
        pdl_flux_c = pdl_flux_c_conv[iii,:]
        out_spec_pdl = interp1d(wave_c, pdl_flux_c,bounds_error=False,fill_value=0.)(wave_unc)
        n_c_out=out_spec_pdl.shape
        error=[]
        for i in range(0,n_unc):
            val=out_spec_pdl[i]
            if np.isnan(val):
                val=0
            model[i][iii]=val*masked[i]
            model_no_mask[i][iii]=val
            if masked[i] > 0:
                error.extend([0.01*abs(e_flux_unc[i])])
            else:
                error.extend([0])
    pdl_C_input=np.zeros(nf)
    pdl_model=np.zeros([n_unc,nf])
    pdl_model_good=np.zeros([n_unc,nf])
    pdl_model_no_mask=np.zeros([n_unc,nf])
    pdl_error=np.zeros(n_unc)
    pdl_masked=masked
    pdl_dust=np.zeros([n_unc,nf])
    for j in range(0, nf):
        for i in range(0, n_unc):
            wave_res=wave_unc[i]/(1+redshift)
            dust_rat=my.A_l(3.1,wave_res)
            dust=10**(-0.4*Av[j]*dust_rat)  
            pdl_dust[i][j]=dust
            val=model[i][j]*dust
            pdl_model[i][j]=val
            val_no_mask=model_no_mask[i][j]*dust
            pdl_model_no_mask[i][j]=val_no_mask
            e_val=error[i]
            val_now=e_flux_unc[i]
            if val_now == 0:
                val_now=1
            pdl_error[i]=1.0/(abs(val_now)**2)
    pdl_flux_masked=flux_masked
#######################################################
# LINEAR GUESS
    j_iter=0
    n_iter=n_mc
    coeffs_iter=np.zeros([nf,n_iter])
    for j_iter in range(0, n_iter):
        pdl_gr=ran.randn(n_unc)
#        pdl_gr->inplace->clip(-1,1);
        pdl_noise=np.sqrt(1/pdl_error)*pdl_gr
        pdl_flux_to_fit=pdl_flux_masked+pdl_noise
        [y_model_now, coeffs] = linfit1d(pdl_flux_to_fit,pdl_model,1/pdl_error)
        pdl_1st_model=y_model_now
        #
        # We remove the models that are negative
        #
        nf_new=0
        nf_neg=0
        for k in range(0, nf):
            C=coeffs[k][0]
            if C > 0:
                nf_new=nf_new+1
            else:
                nf_neg=nf_neg+1
        MOD=[]
        if nf_new > 0:
            while nf_neg > 0:
                pdl_model_new=np.zeros([n_unc,nf_new])
                nf_i=0;
                for k in range(0, nf):
                    C=coeffs[k][0]
                    if C > 0:
                        pdl_model_new[:,nf_i] = pdl_model[:,k]
                        MOD.extend([k])
                        nf_i=nf_i+1
                    else:
                        coeffs[k][0]=0
                [yfit, coeffs_new] = linfit1d(pdl_flux_to_fit,pdl_model_new,1.0/pdl_error)
                pdl_1st_model=yfit
                y_model_now=yfit
                nf_i=0
                nf_neg=0
                nf_new=0
                for k in range(0, nf):
                    C=coeffs[k][0]
                    if C > 0:
                        val=coeffs_new[nf_i][0]
                        nf_i=nf_i+1
                        if val > 0:
                            coeffs[k][0]=val
                            nf_new=nf_new+1
                        else:
                            coeffs[k][0]=0
                            nf_neg=nf_neg+1
                if nf_new == 0:
                    nf_neg=0
            #for k in range(0, nf):
            #    C=coeffs[k][0]
        else:
            nf_new=nf
            #for k in range(0, nf):
            #    C=coeffs[k,0]
        coeffs_iter[:,j_iter] = coeffs[:,0]
    #print coeffs_iter,j_iter
    # End LINEAR GUESS
    #######################################################
    pdl_C_input = np.sum(coeffs_iter,axis=1)
    #print pdl_C_input
    sum_JUNK=np.sum(pdl_C_input)
    pdl_C_input=pdl_C_input/sum_JUNK
    pdl_C_input_zero=np.copy(pdl_C_input)
    coeffs_iter=coeffs_iter/sum_JUNK
    pdl_C_rms=np.copy(pdl_C_input)
    for j in range(0, nf):
        sec=coeffs_iter[j,:]
        stats_sec1=np.std(sec)+np.mean(sec)
        pdl_C_rms[j]=stats_sec1
    min_C=1e12
    max_C=0
    for j in range(0, nf):
        val_C=pdl_C_input_zero[j]
        if val_C > 0:
            if min_C < val_C:
                min_C=val_C
            if max_C > val_C:
                max_C=val_C
    coeffs=np.zeros([nf,3])    
#
# We fit
#
    ini_cat=0
    pdl_C_flast=np.ones(nf)
    fact_q=1
    i_mc=0
    y_model_end=np.zeros(n_unc)
    for j_mc in range(0, n_mc):
        C_left=1
        if i_mc == 0:
            fact_q=0
        if i_mc == 1:
            fact_q=1 
        i_mc=i_mc+1
        pdl_random=ran.rand(nf)
        pdl_grandom=ran.rand(nf)
        pdl_random_J=nf*ran.rand(nf)
        h_nf=1
        jj=0
        y_model_now=np.zeros(n_unc)
        y_model_no_mask_now=np.zeros(n_unc)
        sum_J=0
        pdl_C_now=np.copy(pdl_random)
        if i_mc > 1:
            for j in range(0, nf):
                val_random=pdl_random[j]
                val_grandom=2*pdl_grandom[j]-1;
                C_val_input=pdl_C_input[j]
                C_val_zero=pdl_C_input_zero[j]
                C_val=C_val_zero
                C_rms=pdl_C_rms[j]
                C_val=C_val+2*fact_q*C_rms
                if C_val > 1:
                    C_val=1
                if C_val < 0:
                    C_val=0
                pdl_C_now[j]=C_val
        else:
            pdl_C_now=np.copy(pdl_C_input)
        # We normalize!    
        sum_J=np.sum(pdl_C_now)
        pdl_C_now=pdl_C_now/sum_J
#        print '______A'
#        print pdl_C_now
#        print '______A'
        for j in range(0, nf):
            pdl_model_j=pdl_model[:,j]
            pdl_model_no_mask_j=pdl_model_no_mask[:,j]        
            val=pdl_C_now[j]
            y_model_now=y_model_now+val*pdl_model_j        
            y_model_no_mask_now=y_model_no_mask_now+val*pdl_model_no_mask_j
        dim_now=y_model_now.shape
        j1=int(0.47*n_unc)
        j2=int(0.53*n_unc)
        a_norm=[]
        b_norm=[]
        j_a=0
        band=0
        for j in range(j1, j2):
            if band == 0:
                a_norm.extend([flux_unc[j]])
                b_norm.extend([y_model_now[j]])
                band=1
            else:
                a_norm[j_a]=flux_unc[j]
                b_norm[j_a]=y_model_now[j]
            if a_norm[j_a] > 0 and b_norm[j_a] > 0:
                j_a=j_a+1
                band=0
        med_b=np.median(b_norm)
        if med_b != 0:
            med_norm=np.median(a_norm)/np.median(b_norm)
        else:
            med_norm=1
        MED_NORM.extend([med_norm])
        y_model_now=y_model_now*med_norm
        y_model_no_mask_now=y_model_no_mask_now*med_norm
        ##############################
        # CHISQ VALUE
        chi=0
        chi2=0
        NFREE=0
        have_error=1
        out_spec=np.zeros(n_unc)
        chi_sec=np.zeros(n_unc)
        ran.seed(None)
        pdl_rand_noise=2*ran.rand(n_unc)-1
        for j in range(0, n_unc):
            out_spec[j]=y_model_now[j]
            chi_sec[j]=0
            if flux_unc[j] != 0 and out_spec[j] != 0  and e_flux_unc[j] != 0:
                rnd=0
                chi=chi+masked[j]*((flux_masked[j]+rnd-out_spec[j])**2)/(e_flux_unc[j])**2
                if have_error == 0:
                    chi_sec[j]=masked[j]*((flux_unc[j]+rnd-out_spec[j])**2)/abs(out_spec[j])
                else:
                    chi_sec[j]=masked[j]*((flux_unc[j]+rnd-out_spec[j])**2)/(e_flux_unc[j])**2
                NFREE=NFREE+1
        chi_sq=chi
        if NFREE > 0:
            chi_sq=(chi_sq/(NFREE))**0.5
        out_ps_now="junkMC1";
        title="X = "+str(chi_sq)+" Q="+str(fact_q)
        plot_results(plot,wave_unc,[flux_unc,out_spec,flux_unc-out_spec,y_model_end,e_flux_unc],out_ps_now,title)
        if chi_sq < 1.1*chi_sq_min_now:
            pdl_C_input=pdl_C_now
            fact_q=0.95*fact_q
            j_mc=0
            if fact_q < 0.05 and i_mc > 1:
                j_mc=n_mc        
            if chi_sq < chi_sq_min_now:
                chi_sq_min_now=chi_sq
            coeffs[:,0] = pdl_C_now  
            sum_JUNK=np.sum(pdl_C_now)
            y_model_end=y_model_now
            y_model_no_mask_end=y_model_no_mask_now
            nf_1=nf-1
            #print '______A1'
            #print pdl_C_now[0:nf]
            #print '______A1'
            coeffs_cat[0:nf,ini_cat]=pdl_C_now[0:nf]
            coeffs_cat[nf,ini_cat]=chi_sq
            pdl_model_spec_cat[:,ini_cat]=y_model_no_mask_end
            ini_cat=ini_cat+1
            if ini_cat > n_mc-2:
                j_mc=n_mc
    #
    # We construct the average model
    # and the average coefficients
    #
    #print '____B'
    #print coeffs
    #print '____B'
    #print coeffs_cat
    #print '____C'
    model_spec=np.zeros(n_unc)
    res_spec=np.zeros(n_unc)
    model_spec_min=np.zeros(n_unc)
    out_spec_now=np.zeros(n_unc)
    SUM_W=0
    out_coeffs=np.zeros([nf,1])
    out_coeffs_e=[]
    N_CASES=0
    pdl_model_final=np.zeros(n_unc)
    for J in range(0, ini_cat):
        CHI=coeffs_cat[nf,J]
        if CHI < 1.1*chi_sq_min_now:
            for j in range(0, n_unc):
                out_spec_now[j]=out_spec_now[j]+(pdl_model_spec_cat[j,J])/CHI
            if N_CASES > 0:
                out_coeffs=np.concatenate((out_coeffs,np.zeros([nf, 1])),axis=1)   
            for j in range(0, nf):
                val=coeffs_cat[j,J]
                out_coeffs[j][N_CASES]=val
            N_CASES=N_CASES+1
            SUM_W=SUM_W+1/CHI
    if SUM_W == 0:
        SUM_W=1
        #
        # No better solution found than the 1st one!!!
        #
        for j in range(0, n_unc):
            val=pdl_1st_model[j]
            model_spec[j]=val;
            out_spec[j]=val;
            res_spec[j]=flux_unc[j]-model_spec[j]        
            model_spec_min[j]=model_spec[j]
    else:
        for j in range(0, n_unc):    
            model_spec[j]=out_spec_now[j]/SUM_W
            out_spec[j]=out_spec_now[j]/SUM_W
            res_spec[j]=flux_unc[j]-model_spec[j]        
            model_spec_min[j]=model_spec[j]
    min_coeffs=coeffs
    for j in range(0, nf):
        tmp=np.zeros(N_CASES)
        for J in range(0, N_CASES):
            tmp[J]=out_coeffs[j][J]
        val=np.mean(tmp)
        sigma=np.std(tmp)       
        sigma_MC=pdl_C_rms[j]
        sigma=np.sqrt(sigma**2+sigma_MC**2)
        sum_C=np.sum(pdl_C_input)
        old_val=coeffs[j,0]
        coeffs[j,0]=val
        coeffs[j,1]=sigma
        coeffs[j,2]=old_val
    chi_sq=chi_sq_min_now
    if chi_sq < MIN_CHISQ:
        MIN_CHISQ=chi_sq
        age_min=0
        met_min=0
        Av_min=0
        age_min_mass=0
        met_min_mass=0
        Av_min_mass=0
        norm_C=0
        norm_C_mass=0
        for k in range(0, nf):
            dust=10**(-0.4*Av[k]*dust_rat)
            C=coeffs[k,0]
            norm_C=norm_C+C
            norm_C_mass=norm_C_mass+C*ml[k]
        for k in range(0, nf):
            C=coeffs[k,0]
            #if norm_C > 0:
            coeffs_N[k,0]=C/norm_C
            #if norm_C_mass > 0:
            coeffs_NM[k,0]=C/norm_C_mass
            if norm_C > 0:
                age_min=age_min+C*np.log10(age_mod[k])/norm_C
                met_min=met_min+C*met_mod[k]/norm_C
                Av_min=Av_min+C*Av[k]/norm_C        
                CN=C/norm_C
            C_now=C*med_norm
            if norm_C_mass > 0:
                age_min_mass=age_min_mass+C*np.log10(ml[k]*age_mod[k])/norm_C_mass
                met_min_mass=met_min_mass+C*ml[k]*met_mod[k]/norm_C_mass
                Av_min_mass=Av_min_mass+C*ml[k]*Av[k]/norm_C_mass
        age_min=10**(age_min)
        age_min_mass=10**(age_min_mass)
    #name=unc_file+", "
    scale="1"
    pdl_age_mod=np.array(age_mod)
    pdl_met_mod=np.array(met_mod)
    pdl_ml=np.array(ml)
    pdl_Av=np.array(Av)
    pdl_model_spec_min=np.array(model_spec_min)
    pdl_res=np.array(res_spec)
    out_ps_now="junkMC2"
    title="X="+str(chi_sq)+" Av="+str(Av[0])+" z="+str(redshift)+" sigma="+str(sigma)
    if plot > 1:
        plot_results(plot,wave_unc,[pdl_flux_masked,pdl_model_spec_min,pdl_res],out_ps_now,title)
    return [chi_sq,pdl_age_mod,pdl_met_mod,pdl_ml,pdl_Av,coeffs,coeffs_N,coeffs_NM,pdl_model_spec_min,pdl_res,pdl_C_input_zero,pdl_C_rms] 

def fit_elines_grad_rnd_new_guided(pdl_wave,pdl_flux,pdl_e_flux,n_mod,chi_goal,d_chi_goal,typef,a,ia,a0,a1,link,n_rnd,pdl_masked,deft,SCALE_INI,g_v,g_d):
    a_out=np.zeros([n_mod, 9])
    SCALE=SCALE_INI
    n_mc=n_rnd
    cont=0
    dims=pdl_flux.shape
    if len(dims) == 1:
        NX=dims[0]
    # First Guess
    n_points=np.sum(pdl_masked)
    pdl_model_0=np.zeros(NX)
    pdl_model_cont_0=np.zeros(NX)
    for i in range(0, n_mod):
        pdl_tmp=create_single_model(pdl_wave,i,typef,a)
        pdl_model_0=pdl_model_0+pdl_tmp
    pdl_model_cont_0=create_single_model(pdl_wave,(n_mod-1),typef,a)
    pdl_chi_now=pdl_masked*((pdl_flux-pdl_model_0)**2)/(pdl_e_flux**2)
    chi_sq=np.sum(pdl_chi_now)
    chi_sq=(chi_sq/(n_points-n_mod-1))**0.5
    chi_sq_ini=chi_sq
    a_out=copy_a(n_mod,a)
    a_now=np.zeros([n_mod,9])
    pdl_rnd=ran.randn(n_mod*9*n_mc)
    pdl_rnd_lin=0.8+0.4*ran.rand(n_mod*9*n_mc)
    #
    # 1st we derive the redshift!
    #
    if g_v == 0:
        new_mc=3
    else:
        new_mc=int(n_mc/2.0)
    for ii in range(0, new_mc):
    # We change slighly the parameters
    #
        for i in range(0, n_mod):
            for j in range(0, 9):
                rnd=pdl_rnd[(j+(9*i)+(9*n_mod*ii))]#*ii]
                rnd_lin=pdl_rnd_lin[(j+(9*i)+(9*n_mod*ii))]#*ii]
                if ia[i][j] == 1:
                    if link[i][j] == -1:
                        if typef[i] != "poly1d\n":
                            if j == 3:
                                a_now[i][j]=a0[i][j]+rnd_lin*(a1[i][j]-a0[i][j])#*ii/new_mc
                            else:
                                a_now[i][j]=a_out[i][j]#+$SCALE*$a_out[$i][$j]*$rnd;
                        else:
                            a_now[i][j]=a_out[i][j]+SCALE*a_out[i][j]*rnd
                        if a_now[i][j] < a0[i][j]:
                            a_now[i][j]=a0[i][j]
                        if a_now[i][j] > a1[i][j]:
                            a_now[i][j]=a1[i][j]
                    else:
                        k=link[i][j]-1
                        method=a1[i][j]
                        if method == 0:
                            a_now[i][j]=a_now[k][j]+a0[i][j]
                        else:
                            a_now[i][j]=a_now[k][j]*a0[i][j]
                else:
                    a_now[i][j]=a_out[i][j]
        n_mod_free=0
        for i1 in range(0, n_mod):
            if typef[i1] == "eline\n":
                if ia[i1][1] == 1:
                    if link[i1][1] == -1:
                        pdl_tmp=create_single_model_one(pdl_wave,i1,typef,a_now)
                        if n_mod_free == 0:
                            pdl_model_n=pdl_tmp
                        else:
                            pdl_model_n=np.vstack((pdl_model_n,pdl_tmp))
                        n_mod_free=n_mod_free+1
                    else:
                        pdl_tmp=create_single_model_one(pdl_wave,i1,typef,a_now)
                        ik=link[i1][1]-1
                        t=pdl_model_n[ik,:]
                        pdl_tmp=pdl_tmp*a0[i1][1]        
                        pdl_model_n[ik,:]= t +pdl_tmp           
            if typef[i1] == "poly1d\n":        
                for j in range(0, 9):
                    if ia[i1][j] == 1:
                        pdl_tmp=create_single_model_poly(pdl_wave,j)
                        pdl_model_n=np.vstack((pdl_model_n,pdl_tmp))
                        n_mod_free=n_mod_free+1
        pdl_grandom=ran.randn(NX)
        pdl_flux_fit=pdl_flux-cont+pdl_grandom*pdl_e_flux
        dim_n=pdl_model_n.shape
        if dim_n[1] > 1:
            [pdl_model_0, coeffs_0] = linfit1d(pdl_flux_fit,pdl_model_n)#,{Weights=>1/$pdl_e_flux});#,{Weights=>1/$pdl_e_flux});
        else:
            pdl_model_n[np.where(pdl_model_n == 0)[0]]=1.0
            pdl_rat=pdl_flux_fit/pdl_model_n
            pdl_rat[np.isnan(pdl_rat)]=0
            statst=my.stats(pdl_rat)
            coeffs_0=np.ones(1)
            pdl_model_0=statst[2]*pdl_model_n
            coeffs_0=statst[2]*coeffs_0
        pdl_model_0=pdl_model_0+cont
        n_mod_free=0
        for i1 in range(0, n_mod):
            if typef[i1] == "eline\n":
                if ia[i1][1] == 1:
                    if link[i1][1] == -1:
                        a_now[i1][1]=coeffs_0[n_mod_free]
                        if a_now[i1][1] < a0[i1][1]:
                            a_now[i1][1]=a0[i1][1]
                        if a_now[i1][1] > a1[i1][1]:
                            a_now[i1][1]=a1[i1][1]
                        n_mod_free=n_mod_free+1
                    else:
                        k=link[i1][1]-1
                        method=a1[i1][1]
                        if method == 0:
                            a_now[i1][1]=a_now[k][1]+a0[i1][1]
                        else:
                            a_now[i1][1]=a_now[k][1]*a0[i1][1]
            if typef[i1] == "poly1d\n":
                for j in range(0, 9):
                    if ia[i1][j] == 1:
                        val=coeffs_0[n_mod_free]
                        a_now[i1][j]=val
                        n_mod_free=n_mod_free+1
        pdl_model_0=np.zeros(NX)
        pdl_model_cont_0=np.zeros(NX)
        for i in range(0, n_mod):
            pdl_tmp=create_single_model(pdl_wave,i,typef,a_now)
            pdl_model_0=pdl_model_0+pdl_tmp
        pdl_model_cont_0=create_single_model(pdl_wave,(n_mod-1),typef,a_now)
        pdl_a_now=a_now    
        pdl_chi_now=pdl_masked*((pdl_flux-pdl_model_0)**2)/(pdl_e_flux**2)
        chi_sq=np.sum(pdl_chi_now)
        chi_sq=(chi_sq/(n_points-n_mod-1))**0.5
        if chi_sq <= chi_sq_ini:
            a_out=copy_a(n_mod,a_now)
            chi_sq_ini=chi_sq
            pdl_model=pdl_model_0
            pdl_model_cont=pdl_model_cont_0
    chi_sq_ini=1e12
    a_now=copy_a(n_mod,a_out)
#
# 2nd we derive the sigma
#
    if g_d == 0:
        new_mc=3
    else:
        new_mc=int(n_mc)
    for ii in range(0, new_mc):
        #
        # We change slighly the parameters
        #
        for i in range(0, n_mod):
            for j in range(0, 9):
                rnd=pdl_rnd[(j+(9*i)+(9*n_mod*ii))]#*$ii);
                rnd_lin=pdl_rnd_lin[(j+(9*i)+(9*n_mod*ii))]#*$ii);
                if ia[i][j] == 1:
                    if link[i][j] == -1:
                        if typef[i] != "poly1d\n":
                            if j == 2:
                                a_now[i][j]=a0[i][j]+rnd_lin*(a1[i][j]-a0[i][j])#*$ii/$new_mc
                            else:
                                a_now[i][j]=a_out[i][j]
                        else:
                            a_now[i][j]=a_out[i][j]+SCALE*a_out[i][j]*rnd
                        if a_now[i][j] < a0[i][j]:
                            a_now[i][j]=a0[i][j]
                        if a_now[i][j] > a1[i][j]:
                            a_now[i][j]=a1[i][j]
                    else:
                        k=link[i][j]-1
                        method=a1[i][j]
                        if method == 0:
                            a_now[i][j]=a_now[k][j]+a0[i][j]
                        else:
                            a_now[i][j]=a_now[k][j]*a0[i][j]
                else:
                    a_now[i][j]=a_out[i][j]
        n_mod_free=0
        for i1 in range(0, n_mod):
            if typef[i1] == "eline\n":
                if ia[i1][1] == 1:
                    if link[i1][1] == -1:
                        pdl_tmp=create_single_model_one(pdl_wave,i1,typef,a_now)
                        if n_mod_free == 0:
                            pdl_model_n=pdl_tmp
                        else:
                            pdl_model_n=np.vstack((pdl_model_n,pdl_tmp))
                        n_mod_free=n_mod_free+1
                    else:
                        pdl_tmp=create_single_model_one(pdl_wave,i1,typef,a_now)
                        i2=link[i1][1]-1
                        t=pdl_model_n[i2,:]
                        pdl_tmp=pdl_tmp*a0[i1][1]       
                        pdl_model_n[i2,:]= t +pdl_tmp 
            if typef[i1] == "poly1d\n":        
                for j in range(0, 9):
                    if ia[i1][j] == 1:
                        pdl_tmp=create_single_model_poly(pdl_wave,j)
                        pdl_model_n=np.vstack((pdl_model_n,pdl_tmp))
                        n_mod_free=n_mod_free+1
        pdl_grandom=ran.randn(NX)
        pdl_flux_fit=pdl_flux-cont+pdl_grandom*pdl_e_flux
        dim_n=pdl_model_n.shape
        if dim_n[1] > 1:
            [pdl_model_0, coeffs_0] = linfit1d(pdl_flux_fit,pdl_model_n)
        else:
            pdl_model_n[np.where(pdl_model_n == 0)[0]]=1.
            pdl_rat=pdl_flux_fit/pdl_model_n
            pdl_rat[np.isnan(pdl_rat)]=0
            statst=my.stats(pdl_rat)
            coeffs_0=np.ones(1)
            pdl_model_0=statst[2]*pdl_model_n
            coeffs_0=statst[2]*coeffs_0
        pdl_model_0=pdl_model_0+cont
        n_mod_free=0
        for i1 in range(0, n_mod):
            if typef[i1] == "eline\n":
                if ia[i1][1] == 1:
                    if link[i1][1] == -1:
                        a_now[i1][1]=coeffs_0[n_mod_free]
                        n_mod_free=n_mod_free+1
                    else:
                        k=link[i1][1]-1
                        method=a1[i1][1]
                        if method == 0:
                            a_now[i1][1]=a_now[k][1]+a0[i1][1]
                        else:
                            a_now[i1][1]=a_now[k][1]*a0[i1][1]
            if typef[i1] == "poly1d\n":
                for j in range(0, 9):
                    if ia[i1][j] == 1:
                        val=coeffs_0[n_mod_free]
                        a_now[i1][j]=val
                        n_mod_free=n_mod_free+1
        pdl_model_cont_0=create_single_model(pdl_wave,(n_mod-1),typef,a_now)
        pdl_a_now=a_now
        pdl_chi_now=pdl_masked*((pdl_flux-pdl_model_0[:,0])**2.0)/(pdl_e_flux**2)
        chi_sq=np.sum(pdl_chi_now)
        chi_sq=(chi_sq/(n_points-n_mod-1))**0.5
        if chi_sq <= chi_sq_ini:
            a_out=copy_a(n_mod,a_now)
            chi_sq_ini=chi_sq
            pdl_model=pdl_model_0
            pdl_model_cont=pdl_model_cont_0
    chi_sq_ini=1e12
    a_now=copy_a(n_mod,a_out)
    # We start the fitting loop!
    ii=0
    SCALE_IN=SCALE
    new_mc=int(n_mc/3)
    while ii < new_mc:
    #
    # We change slighly the parameters
    #
        for i in range(0, n_mod):
            for j in range(0, 9):
                rnd=pdl_rnd[(j+(9*i)+(9*n_mod*ii))]#*$ii)
                rnd_lin=pdl_rnd_lin[(j+(9*i)+(9*n_mod*ii))]#*$ii);
                if ia[i][j] == 1:
                    if link[i][j] == -1:
                        A1=a1[i][j]
                        A0=a0[i][j]
                        if a1[i][j] > 1.3*a_out[i][j]:
                            A1=1.3*a_out[i][j]
                        if a0[i][j] < 0.7*a_out[i][j]:
                            A0=0.7*a_out[i][j]
                        if typef[i] == "eline\n":
                            if j == 3 or j == 2:
                                a_now[i][j]=a_out[i][j]+SCALE*rnd*(A1-A0)/(5*new_mc)
                            else:
                                a_now[i][j]=a_out[i][j]+SCALE*a_out[i][j]*rnd
                            if a_now[i][j] < A0:
                                a_now[i][j]=A0
                            if a_now[i][j] > A1:
                                a_now[i][j]=A1
                        else:
                            # No variation!
                            a_now[i][j]=a_out[i][j]#+$SCALE*$a_out[$i][$j]*$rnd*0.0001;
                    else:
                        k=link[i][j]-1
                        method=a1[i][j]
                        if method == 0:
                            a_now[i][j]=a_now[k][j]+a0[i][j]
                        else:
                            a_now[i][j]=a_now[k][j]*a0[i][j]
                else:
                    a_now[i][j]=a_out[i][j]
        for i1 in range(0, n_mod):
            if typef[i1] == "eline\n":
                if ia[i1][1] == 1:
                    if link[i1][1] == -1:
                        if a_now[i1][1] < a0[i1][1]:
                            a_now[i1][1]=a0[i1][1]
                        if a_now[i1][1] > a1[i1][1]:
                            a_now[i1][1]=a1[i1][1]
        # MC test
        pdl_model_0=np.zeros(NX)
        pdl_model_cont_0=np.zeros(NX)
        for i in range(0, n_mod):
            pdl_tmp=create_single_model(pdl_wave,i,typef,a_now)
            pdl_model_0=pdl_model_0+pdl_tmp
        pdl_model_cont_0=create_single_model(pdl_wave,(n_mod-1),typef,a_now)
        ########################
        # 2014.11.11 LAST CHANGE
        #
        pdl_grandom=ran.randn(NX)
        pdl_flux_fit=pdl_flux-cont+pdl_grandom*pdl_e_flux
        pdl_chi_now=pdl_masked*((pdl_flux_fit-pdl_model_0)**2)/(pdl_e_flux**2)
        chi_sq=np.sum(pdl_chi_now)
        chi_sq=(chi_sq/(n_points-n_mod-1))**0.5
        if chi_sq < chi_sq_ini:
            a_out=copy_a(n_mod,a_now)
            chi_sq_ini=chi_sq
            pdl_model=pdl_model_0
            pdl_model_cont=pdl_model_cont_0
            SCALE=SCALE*0.99
            if SCALE < 0.1*SCALE_IN:
                SCALE=SCALE_IN*0.1
        else:
            SCALE=SCALE_IN
            if abs(chi_sq-chi_sq_ini) < d_chi_goal or chi_sq_ini < chi_sq_goal:
                ii=n_mc
        ii=ii+1
    chi_sq_now=chi_sq_ini
    return [a_out,chi_sq_now,pdl_model]

