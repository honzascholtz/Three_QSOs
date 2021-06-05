#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 17 14:29:47 2020

@author: jansen
"""

#importing modules
import numpy as np
import matplotlib.pyplot as plt; plt.ioff()

from astropy.io import fits as pyfits
from astropy import wcs
from astropy.table import Table, join, vstack
from matplotlib.backends.backend_pdf import PdfPages
import pickle
from scipy.optimize import curve_fit

import Graph_setup as gst 

from scipy.misc import imresize

nan= float('nan')

pi= np.pi
e= np.e

plt.close('all')
c= 3.*10**8
h= 6.62*10**-34
k= 1.38*10**-23

Ken98= (4.5*10**-44)
Conversion2Chabrier=1.7 # Also Madau
Calzetti12= 2.8*10**-44
arrow = u'$\u2193$' 


PATH='/Users/jansen/Google Drive/Astro/'
fsz = gst.graph_format()

Sample = Table.read(PATH+'Four_Quasars/Four_Quasars.fits')
New = Sample

itere= np.array([0,1,2])

jj=2

import Plotting_tools as emplot
import Fitting_tools as emfit
import IFU_tools_QSO as IFU

#Plots =  PdfPages(PATH+'Four_quasars/Graphs/Spectra_model_plus_outflow_1G_ll.pdf')

for i in itere:
    ID = Sample['ID'][i].strip()
    
    z = Sample['z'][i]
    
# =============================================================================
#  IMAGES
# ==========================================================================

    if i==3:
        
        Ra_opt = Sample['RA'][i]
        Dec_opt = Sample['DEC'][i]
    
    else:
        Gaia = Table.read(PATH+'Four_Quasars/Catalogues/'+ID+'_Gaia.tbl' , format='ipac')
    
    
        Ra_opt = float(Gaia['ra'])
        Dec_opt = float(Gaia['dec'])
        
    
    
    if ID=='XID_2028':
        storage = emplot.load_obj('KMOS_SIN/Results_storage/Props/lid_1565_H')
    else:
        storage = emplot.load_obj('KMOS_SIN/Results_storage/Props/'+ID+'_H')
        
    storage['Cat_entry'] = New[i]
    storage['X-ray ID'] = ID
    
    

        

# =============================================================================
#         Spectra
# =============================================================================  
    if ID=='HB89':
        IDp='HB8903'
    if ID=='2QZJ':
        IDp='2QZJ00'       
    if ID=='LBQS':
        IDp= 'LBQS01'
    ################
    # Spectrum - Halpha
    ################
    
    #ax1 = axes[ii,jj]
    f, ax1 = plt.subplots(1)
    #ax1 = f.add_axes((0.0625, 0.745-jj*0.325 ,0.425 ,0.25))
    #axres = f.add_axes((0.0625, 0.695-jj*0.325 ,0.425 ,0.05))
    #plt.setp(ax1.get_yticklabels()[:], visible=False)  
    #axres.tick_params(direction='in')
    
    #plt.setp(ax1.get_xticklabels()[:], visible=False) 
    ax1.tick_params(direction='in')
    #if i==0:
    #    xwaaaa=0
    #else:
    #    plt.setp(ax1.get_xticklabels()[:], visible=False) 
    
    
    Data = np.loadtxt(PATH+'Four_Quasars/Results_storage/Spectrums/'+ID+'_H_inn.txt')
    
    wave = Data[0,:]
    flux = np.ma.array(data=Data[1,:], mask=Data[3,:])
    error = Data[2,:]    

    
    #out ,chi2 = emfit.fitting_Halpha_mul_bkp_only(wave,flux,error,z, broad=1)
    
    if i!=2:
        
        out_old, chi2_old = emfit.fitting_Halpha_mul_LBQS(wave,flux,error,z, broad=1)
        
        out, chi2 = emfit.fitting_Halpha_mul_outflow_1G(wave,flux,error,z, broad=1)
        #out, chi2 =    emfit.fitting_Halpha_mul_HB89_narrow(wave,flux,error,z, broad=1)
        
        out_dd, chi2_dd = emfit.fitting_Halpha_mul_LBQS_testing(wave,flux,error,z, broad=1)
        
        print ('BIC between shifted and vel', out.bic-out_dd.bic)
    
    else:
        out, chi2 = emfit.fitting_Halpha_mul_HB89_outflow_1G(wave,flux,error,z, broad=1)
        
        out_old, chi2old =   emfit.fitting_Halpha_mul_testing(wave,flux,error,z, broad=1)
    
    
    
    print ('New BIC ', out.bic)
    print ('Old BIC ' , out_old.bic)
    print (out.bic-out_old.bic)
    
    
    
    print (( out.params['Hao_center'].value-out.params['Han_center'].value)/out.params['Han_center'].value*3e5)
    print (( out.params['Han_fwhm'].value)/out.params['Hao_center'].value*3e5)

    
    z =  (out.params['Han_center'].value*1e4/6562.)-1
    
    
    emplot.plotting_Halpha(wave, flux, ax1, out, 'mul',z, title=0)
    
    ax1.set_ylim(-0.1*max(out.eval(x=wave)), 1.2*max(out.eval(x=wave)))
    
    
    ax1.fill_between((6710, 6736), y1 =-0.1*max(out.eval(x=wave)), y2 = 1.2*max(out.eval(x=wave)) , color='grey', alpha=0.4)
    
    ax1.text(6737, max(out.eval(x=wave)), '[SII]', color='firebrick')
    
    ax1.set_xlim(6300,6900)
    
    ax1.text(6320, 1.1*max(out.eval(x=wave)), IDp+r' - H$\alpha$', fontsize=15)
    ax1.text(6320, 1.0*max(out.eval(x=wave)), r'$\delta$BIC= %.1f' %(out.bic-out_old.bic), fontsize=15)
    
    
    ax1.plot(wave/(1+z)*1e4, out.eval_components(x=wave)['X_'], 'k--' )
    
    try:
        
        ax1.plot(wave/(1+z)*1e4, out.eval_components(x=wave)['Hao_'], 'r--')
        
        try:
            ax1.plot(wave/(1+z)*1e4, out.eval_components(x=wave)['Nro_'], 'r--')
            ax1.plot(wave/(1+z)*1e4, out.eval_components(x=wave)['Nbo_'], 'r--')
            
        except:
            print('No out [NII]')
    except:
        print('No out')
   
    ax1.set_ylabel(r'x$10^{-17}$ (erg/s/cm$^2$/ang)', fontsize=15)
    
    if jj==2:
        ax1.set_xlabel('Rest wavelength (ang)', fontsize=15)
        
    
    wv_rest = wave/(1+z)*1e4
      
    flux_pl = flux.data[np.invert(flux.mask)]
    wv_rst_sc= wv_rest[np.invert(flux.mask)]
    
    fit_loc = np.where((wv_rst_sc>6200)&(wv_rst_sc<7000))[0]
    
    
    Plots.savefig(f)
    
    
    g, axres = plt.subplots(1)
    axres.plot(wv_rst_sc[fit_loc], flux_pl[fit_loc]-out.eval(x=wv_rst_sc[fit_loc]/1e4*(1+z)), drawstyle='steps-mid')
    axres.axhline(0, linestyle='dashed', color='black', alpha=0.5)
    
    axres.fill_between((6710, 6736), y1 =-5,  y2 = 5 , color='grey', alpha=0.4)
    
    axres.set_ylim(-1.2,1.2)
    
    axres.set_xlim(6400,6690)
    axres.set_xlabel('Rest wavelength (ang)', fontsize=15)
    axres.set_ylabel(r'x$10^{-17}$ (erg/s/cm$^2$/ang)', fontsize=15)   
    
    #Plots.savefig(g)
    
#Plots.close()

    
plt.show()