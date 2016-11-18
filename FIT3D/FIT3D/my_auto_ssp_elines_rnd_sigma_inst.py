#!/usr/bin/python

import sys
import numpy as np
from scipy import linalg as li
from numpy import float_
from numpy import absolute as abs
from numpy import random as ran
import matplotlib
from scipy.signal.signaltools import convolve2d
from scipy.ndimage.filters import convolve1d
from scipy.interpolate.interpolate import interp1d
import my as my
import my_auto_ssp_elines_rnd as ssp

def fit_ssp_lin_no_zero(redshift,sigma,Av_NOW,crval,cdelt,crpix,nf,n_c,pdl_flux_c_ini,hdr_c_ini,wave_unc,masked,e_flux_unc,flux_unc,n_mc,chi_sq_min_now,MIN_CHISQ,plot,sigma_inst):
    Av=Av_NOW
    iter_max=5
    last_chi=1e12
    n_unc=len(flux_unc)#+1
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
    #rsigma=sigma/dpix_c

    name=[]
    age_mod=[]
    met_mod=[]
    ml=[]
    flux_c=[]
    pdl_flux_c_conv=[]
    model=np.zeros([n_unc,nf])
    model_no_mask=np.zeros([n_unc,nf])  
    y_model_end=np.zeros([n_unc,nf])
    y_model_no_mask_end=np.zeros([n_unc,nf]) 
    MED_NORM=0
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
        pdl_flux_c_ini_now=pdl_flux_c_ini[iii,:]
        out_spec_pdl=shift_convolve_pdl(np.array(wave_unc),np.array(wave_c),pdl_flux_c_ini_now,0,sigma_inst,sigma)
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
    [y_model_now, coeffs] = ssp.linfit1d(pdl_flux_masked,pdl_model,1.0/pdl_error)
    y_model_now=y_model_now[:,0]
    #print Av,"HOLA"
    #print coeffs
    #sigma_mean=np.sqrt(sigma_inst**2.0+(5000.0*(sigma/300000.0))**2.0)
    #smooth_ratio=ssp.smooth_ratio(np.array(flux_unc),y_model_now,int(sigma_mean))
    #y_model_now=y_model_now*smooth_ratio
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
            [yfit, coeffs_new] = ssp.linfit1d(pdl_flux_masked,pdl_model_new,1.0/pdl_error)
            y_model_now=yfit[:,0]
            #sigma_mean=np.sqrt(sigma_inst**2.0+(5000.0*(sigma/300000.0))**2.0)
            #print sigma_mean,sigma_inst,sigma,y_model_now.shape
            #smooth_ratio=ssp.smooth_ratio(np.array(flux_unc),yfit[:,0],int(sigma_mean))
            #yfit=yfit[:,0]*smooth_ratio
            #y_model_now=yfit
            #print y_model_now.shape
            
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
#    y_model_no_mask_end=y_model_no_mask_now
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
        #sigma_mean=np.sqrt(sigma_inst**2.0+(5000.0*(sigma/300000.0))**2.0)
        #smooth_ratio=ssp.smooth_ratio(np.array(flux_unc),pdl_model_spec_min,int(sigma_mean))
        #pdl_model_spec_min=pdl_model_spec_min*smooth_ratio   
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
    #print title,NFREE
    if plot > 0:
        ssp.plot_results(plot,pdl_wave_unc,[pdl_flux_masked,pdl_model_spec_min,pdl_res],out_ps_now,title)
    return [chi_sq,pdl_age_mod,pdl_met_mod,pdl_ml,pdl_Av,coeffs,coeffs_N,coeffs_NM,pdl_model_spec_min,pdl_res]

def fit_ssp_lin_no_zero_no_cont(redshift,sigma,Av_NOW,crval,cdelt,crpix,nf,n_c,pdl_flux_c_ini,hdr_c_ini,wave_unc,masked,e_flux_unc,flux_unc,n_mc,chi_sq_min_now,MIN_CHISQ,plot,sigma_inst):
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
    #rsigma=sigma/dpix_c
    name=[]
    age_mod=[]
    met_mod=[]
    ml=[]
    flux_c=[]
    model=np.zeros([n_unc,nf])
    model_no_mask=np.zeros([n_unc,nf])
   # y_model_end=
   # y_model_no_mask_end=
   # MED_NORM=    
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
        pdl_flux_c_ini_now=pdl_flux_c_ini[iii,:]
        out_spec_pdl=shift_convolve_pdl(np.array(wave_unc),np.array(wave_c),pdl_flux_c_ini_now,0,sigma_inst,sigma)
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
            pdl_error[i]=1.0/(abs(val_now)**2.0)
#
# We fit
#
    pdl_flux_masked=flux_masked
    ini_cat=0
    pdl_C_flast=np.ones(nf)
# Just a linear FIT, without restrictions!
    [y_model_now, coeffs] = ssp.linfit1d(pdl_flux_masked,pdl_model,1.0/pdl_error)
    sigma_mean=np.sqrt(sigma_inst**2.0+(5000.0*(sigma/300000.0))**2.0)
    smooth_rati=ssp.smooth_ratio(flux_unc,y_model_now[:,0],int(sigma_mean))
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
            [yfit, coeffs_new] = ssp.linfit1d(pdl_flux_masked,pdl_model_new,1.0/pdl_error)
            sigma_mean=np.sqrt(sigma_inst**2.0+(5000.0*(sigma/300000.0))**2.0)
            smooth_rati=ssp.smooth_ratio(flux_unc,yfit[:,0],int(sigma_mean))
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
    #y_model_no_mask_end=y_model_no_mask_now
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
    sigma_mean=np.sqrt(sigma_inst**2.0+(5000.0*(sigma/300000.0))**2.0)  
    smooth_rati=ssp.smooth_ratio(flux_unc,pdl_model_spec_min,int(sigma_mean))
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
        ssp.plot_results(plot,pdl_wave_unc,[pdl_flux_masked,y_model_end,pdl_res],out_ps_now,title)
        #sys.exit()
    return [chi_sq,pdl_age_mod,pdl_met_mod,pdl_ml,pdl_Av,coeffs,coeffs_N,coeffs_NM,pdl_model_spec_min,pdl_res]

def fit_ssp_lin_MC(redshift,sigma,Av_NOW,crval,cdelt,crpix,nf,n_c,pdl_flux_c_ini,hdr_c_ini,wave_unc,masked,e_flux_unc,flux_unc,n_mc,chi_sq_min_now,MIN_CHISQ,plot,sigma_inst):
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
    #rsigma=sigma/dpix_c
    name=[]
    age_mod=[]
    met_mod=[]
    ml=[]
    flux_c=[]
    model=np.zeros([n_unc,nf])
    model_no_mask=np.zeros([n_unc,nf])
    y_model_end=np.zeros([n_unc,nf])
    y_model_no_mask_end=np.zeros([n_unc,nf])
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
        pdl_flux_c_ini_now=pdl_flux_c_ini[iii,:]
        out_spec_pdl=shift_convolve_pdl(np.array(wave_unc),np.array(wave_c),pdl_flux_c_ini_now,0,sigma_inst,sigma)        
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
        [y_model_now, coeffs] = ssp.linfit1d(pdl_flux_to_fit,pdl_model,1.0/pdl_error)
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
                [yfit, coeffs_new] = ssp.linfit1d(pdl_flux_to_fit,pdl_model_new,1.0/pdl_error)
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
        ssp.plot_results(plot,wave_unc,[flux_unc,out_spec,flux_unc-out_spec,y_model_end,e_flux_unc],out_ps_now,title)
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
        ssp.plot_results(plot,wave_unc,[pdl_flux_masked,pdl_model_spec_min,pdl_res],out_ps_now,title)
    return [chi_sq,pdl_age_mod,pdl_met_mod,pdl_ml,pdl_Av,coeffs,coeffs_N,coeffs_NM,pdl_model_spec_min,pdl_res,pdl_C_input_zero,pdl_C_rms] 

def shift_convolve_pdl(pdl_wave_out,pdl_wave_in,pdl_flux_in,redshift,sigma_inst,sigma_km_h):
    #sigma_km_h=150.0
    sigma_inst=sigma_inst*(1+redshift)
    sigma_km_h=sigma_km_h/(1+redshift)
    sub_pix=1.0#+ran.rand()
    wmin=np.amin(np.array(pdl_wave_in))
    wmax=np.amax(np.array(pdl_wave_in))
    dpix_ini=pdl_wave_in[1]-pdl_wave_in[0]
    dpix_out=pdl_wave_out[1]-pdl_wave_out[0]
    dpix=dpix_ini/sub_pix
    n_sub=(wmax-wmin)/dpix
    pdl_wave=wmin+((wmax-wmin)/n_sub)*np.arange(0,n_sub)
    pdl_flux=interp1d(pdl_wave_in,pdl_flux_in,bounds_error=False,fill_value=0.)(pdl_wave)
    rsigma=sigma_inst/dpix
    #rsigma=sigma_km_h/dpix
    e=np.exp(1)
    # Convolved with Sigma_inst
    #pdl_flux_conv_inst
    #print dpix_ini
    #print dpix
    #print sigma_inst
    #print rsigma
    if sigma_inst > 0.1:
        box=int((n_sub-1)/2)+1
        box=int(3.0*rsigma*2.0)
        #print box
        if box < 3:
            box=3
        kernel=np.zeros(2*box+1)
        norm=0
        for j in range(0, 2*box+1):
            gaus=np.exp(-0.5*(((j-box)/rsigma)**2.0))    
            kernel[j]=gaus
            norm=norm+gaus
        kernel=kernel/norm
        #print kernel
        #print pdl_flux.shape
        pdl_flux_conv_inst = convolve1d(pdl_flux,kernel)#,mode='same')
    else:
        pdl_flux_conv_inst=pdl_flux
    l_pdl_wave=np.log(pdl_wave)
    l_pdl_wave0=l_pdl_wave[0]
    l_pdl_wave1=l_pdl_wave[len(pdl_wave)-1]
    # Take care: Factor 5!!!!
    # Interpolation Factor
    f_fine=1.0#5
    d_l_pdl_wave=(l_pdl_wave[1]-l_pdl_wave[0])/f_fine
    n_l_pdl_wave=int((l_pdl_wave1-l_pdl_wave0)/(d_l_pdl_wave))
    pdl_ones=np.arange(0,n_l_pdl_wave)
    l_pdl_wave_l=l_pdl_wave0+d_l_pdl_wave*pdl_ones
    pdl_wave_l=e**(l_pdl_wave_l)/(1+redshift)
    pdl_flux_l_wave=interp1d(pdl_wave,pdl_flux_conv_inst,bounds_error=False,fill_value=0.)(pdl_wave_l)
    ######################################################
    # Convolved with Sigma_km_h
    #
    if sigma_km_h == 0:
        sigma_km_h=1
    rsigma_km_h=((sigma_km_h/300000.0)/(d_l_pdl_wave))
    #print rsigma_km_h,"HOLA",sigma_km_h,d_l_pdl_wave,l_pdl_wave[1],l_pdl_wave[0]
    box=int(5.0*(500.0/300000.0)/(d_l_pdl_wave))
    if box < 3:
        box=3
    kernel=np.zeros(2*box+1)
    norm=0
    for j in range(2*box+1):
        gaus=np.exp(-0.5*(((j-box)/rsigma_km_h)**2))    
        kernel[j]=gaus
        norm=norm+gaus
    kernel=kernel/norm
    #
    # Check THIS!!! not wave_l???
    pdl_flux_conv_km_h =  convolve1d(pdl_flux_l_wave,kernel)#,mode='same') pdl_flux_l_wave#
    #
    #
    ###################################################
    pdl_wave_l_redshift=pdl_wave_l
    pdl_flux_l_wave_km_h=interp1d(pdl_wave_l_redshift,pdl_flux_conv_km_h,bounds_error=False,fill_value=0.)(pdl_wave_out)
    return pdl_flux_l_wave_km_h