#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon May 20 12:58:41 2019

@author: jansen
"""

#importing modules
import numpy as np
import matplotlib.pyplot as plt

from astropy.io import fits as pyfits
from astropy import wcs
from astropy.table import Table, join, vstack
from matplotlib.backends.backend_pdf import PdfPages
import pickle
from scipy.optimize import curve_fit

import Graph_setup as gst 

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



import Graph_setup as gst 
import Tools_IFU as IFU
import Tools_plotting as emplot
import Tools_fitting as emfit
import Tools_path as ph

fsz = gst.graph_format()


plot_it = 0
#0,4

Sample = Table.read(ph.MyPATH+'Four_Quasars.fits')

itere=np.array([1])
for i in itere:    
    ID = Sample['ID'][i].strip()
    
    z= Sample['z'][i]
    Hal_band='K_sin'
        
    H_file = Sample['Halpha cube'][i].strip() 
        
     
    storage_H = IFU.load_cube(ph.MyPATH+'SINFONI/'+H_file+'.fits', z, ID, 'Sinfoni')
    
    storage_H = IFU.add_res(storage_H, Sample[i])
    
    # Masking Emission Region
    storage_H = IFU.mask_emission(storage_H, z)

    # Very quick Sky masking
    storage_H = IFU.mask_sky(storage_H, 1.5)
    
    # Unmasking emission region
    storage_H = IFU.unmask_em(storage_H,z)

    
    # Collapsing the cube to image
    storage_H = IFU.collapse_white(storage_H, plot_it)
        
        
    storage_H = IFU.find_center(storage_H,'Median_stack_white', plot_it)
    

       
        
    storage_H = IFU.choose_pixels(storage_H, plot_it)
        
    storage_H = IFU.astrometry_correction_GAIA(storage_H)
        
        
        
    storage_H = IFU.stack_sky(storage_H, Hal_band, plot_it, expand=0)
    
    storage_H = IFU.D1_spectra_collapse(storage_H,z, 'H', plot_it, addsave='_tot')
        
    storage_H= IFU.fitting_collapse_Halpha(storage_H,z, 1, broad=1)
    
    
    rads = np.linspace(0.1,2,20)
    
    flux_save = np.zeros((4, len(rads)))
    
    #flux_save = np.loadtxt(ph.MyPATH+'KMOS_SIN/Results_storage/Growth/'+ID+'Growth_Hal.txt')
    flux_mes =  np.zeros_like(rads) #flux_save[1,:]
    
    flux_BLR = np.zeros_like(rads)
    
    #flux_mes_old = flux_mes.copy()
    errors = np.zeros(len(rads))
    Hal_plots = PdfPages(ph.MyPATH+'Graphs/'+ID+'_Halpha_ap.pdf')
    
    storage = storage_H
    
    arc = range(len(rads))
    #arc = np.array([9])
    for j in arc:
    
        storage = storage_H
        flx = storage['flux'].data
        center = storage['Median_stack_white_Center_data'][1:3].copy()
        shapes = storage['dim']
    
        mask_catch = storage['flux'].mask.copy()
        mask_catch[:,:,:] = True
        header  = storage['header']
        
        try:
            arc = np.round(1./(header['CD2_2']*3600))
        except:
            arc = np.round(1./(header['CDELT2']*3600))
            
    
        rad_sel = arc* rads[j]
    
        print ('radius ', rad_sel, ' pixels, ', rads[j], ' arcseconds')
        
        # This choose spaxel within certain radius. Then sets it to False since we dont mask those pixels
        for ix in range(shapes[0]):
            for iy in range(shapes[1]):
                dist = np.sqrt((ix- center[1])**2+ (iy- center[0])**2)
                if dist< rad_sel:
                    mask_catch[:,ix,iy] = False
    
        mask_sky_1D = storage['sky_clipped_1D'].copy()
    
        flux = np.ma.array(data=storage['flux'].data, mask= mask_catch) 
        
        D1_spectra = np.ma.sum(flux, axis=(1,2))
        D1_spectra = np.ma.array(data = D1_spectra.data, mask=mask_sky_1D)
        
        
        wave = storage['obs_wave'].copy()
        error = IFU.STD_calc(wave/(1+z)*1e4,D1_spectra, 'H')* np.ones_like(D1_spectra)
        
        try:
            D_out = storage['1D_fit_Halpha_mul']
            loc = D_out.params['Han_center'].value 
            
            #loc_w = D_out.params['Haw_center'].value 
            #Hal_cm = 6562.*(1+z)/1e4
            #wid_w = (D_out.params['Haw_fwhm'].value)/Hal_cm*2.9979e5
        
            #dec = np.array([wid_w,loc_w])
            
            out, chi2 = emfit.fitting_Halpha_mul_bkp(wave,D1_spectra,error,z,broad=1,decompose=D_out )        
    
        
            
            
            f, (ax1,ax2) = plt.subplots(2)      
            ax1.imshow(np.ma.array(data=storage['Median_stack_white'], mask=mask_catch[0,:,:]), origin='low')     
            emplot.plotting_Halpha(wave, D1_spectra, ax2, out, 'mul',z)
            
            
            
            fl, error = IFU.flux_measure_ind(out,wave, 'H', use='BPT', error=1)   
            fl_blr, error_blr = IFU.flux_measure_ind(out,wave, 'H', use='broad', error=1)   
            
            flux_mes[j] = fl
            errors[j] = error
            
            flux_BLR[j] = fl_blr
            print ('flux and error ', fl, error)
            ax1.set_title('Rad '+str(rads[j])+ ', flx '+str(fl) +' ' + str(fl_blr))
    
            Hal_plots.savefig()
        except:
            print ('?')
    
    Hal_plots.close()
    g,axflx = plt.subplots(1)   
    axflx.errorbar(rads, flux_mes/max(flux_mes),yerr= errors/max(flux_mes), label='Halpha', color='red')
    axflx.errorbar(rads, flux_BLR/max(flux_BLR),yerr= 0.1, label='BLR', color='blue')
    
    
    flux_save = np.zeros((5, len(rads)))
    flux_save[0,:] = rads
    flux_save[1,:] = flux_mes
    #flux_save[2,:] = PSF
    flux_save[3,:] = errors
    flux_save[4,:] = flux_BLR
    
    plt.show()
    
    #np.savetxt(ph.MyPATH+'Results_storage/Growth/'+ID+'Growth_Hal.txt', flux_save)
    
    plt.xlabel('Radius (arcseconds)')
    plt.ylabel('Normalized flux')
    
    #plt.close('all')

plt.show()