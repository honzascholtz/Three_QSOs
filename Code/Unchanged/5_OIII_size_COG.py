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


PATH='/Users/jansen/Google Drive/Astro/'
fsz = gst.graph_format()



binning = 'Nearest'



import Graph_setup as gst 
import IFU_tools_QSO as IFU
import Plotting_tools as emplot
import Fitting_tools as emfit


plot_it = 0
#0,4

Sample = Table.read(PATH+'Four_Quasars/Four_Quasars.fits')
Sample['Halpha cube'][0] = '2QZ0028-28'


OIII_size={'HB89': np.array(['0.2', '1'], dtype=str)}
OIII_size['LBQS'] = np.array(['0.2', '1'], dtype=str)
OIII_size['2QZJ'] = np.array(['0.2', '1'], dtype=str)

itere=np.array([2])
for i in itere:    
    ID = Sample['ID'][i].strip()
    H_file=ID+'_H_fl'
    z= Sample['z'][i]
    OIII_band='Hsin'
        
    if ID=='2QZJ':
        z= 2.4064
        
    storage_O = storage_H = IFU.load_cube(PATH+'Four_Quasars/SINFONI/'+H_file+'.fits', z, ID, 'Sinfoni')

    storage_O = IFU.add_res(storage_O, Sample[i])

    # Masking Emission Region
    storage_O = IFU.mask_emission(storage_O, z)

    # Very quick Sky masking
    storage_O = IFU.mask_sky(storage_O, 1.5)

    # Unmasking emission region
    storage_O = IFU.unmask_em(storage_O,z)


    # Collapsing the cube to image
    storage_O= IFU.collapse_white(storage_O, 0)

    storage_O = IFU.find_center(storage_O,'Median_stack_white', 0)

    size = OIII_size[ID]

    rds = float(size[0])
    fl = size[1]

    print ('Radius and extraction region ', rds, fl)
    storage_O = IFU.choose_pixels(storage_O, plot_it, rad= rds , flg = fl)

    storage_O = IFU.astrometry_correction_GAIA(storage_O)

    storage_O = IFU.stack_sky(storage_O, OIII_band, 0, expand=0)
    
    storage_O = IFU.D1_spectra_collapse(storage_O,z, 'OIII', 0, addsave='_tot')

    storage_O= IFU.fitting_collapse_OIII(storage_O,z, 1)
    
    
    rads = np.linspace(0.1,2,20)
    
    flux_save = np.zeros((4, len(rads)))
    
    #flux_save = np.loadtxt(PATH+'KMOS_SIN/Results_storage/Growth/'+ID+'Growth_Hal.txt')
    flux_mes =  np.zeros_like(rads) #flux_save[1,:]
    
    flux_BLR = np.zeros_like(rads)
    
    #flux_mes_old = flux_mes.copy()
    errors = np.zeros(len(rads))
    Hal_plots = PdfPages(PATH+'Four_Quasars/Graphs/'+ID+'_OIII_ap_COG.pdf')
    
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
        error = IFU.STD_calc(wave/(1+z)*1e4,D1_spectra, 'OIII')* np.ones_like(D1_spectra)
        
        if 1==1:
            outo = storage['1D_fit_OIII_mul']
            
            out = emfit.fitting_OIII_Hbeta_qso_mul(wave,D1_spectra, error,z, decompose=outo)
    
        
            
            
            f, (ax1,ax2) = plt.subplots(2)      
            ax1.imshow(np.ma.array(data=storage['Median_stack_white'], mask=mask_catch[0,:,:]), origin='low')     

            emplot.plotting_OIII_Hbeta(wave, D1_spectra, ax2, out, 'mul',z, title=0)
            
            
            fl, error = IFU.flux_measure_ind(out,wave, 'OIII', use='tot', error=1)   
            fl_blr, error_blr = IFU.flux_measure_ind(out,wave, 'Hb', use='bro', error=1)   
            
            flux_mes[j] = fl
            errors[j] = error
            
            flux_BLR[j] = fl_blr
            print ('flux and error ', fl, error)
            ax1.set_title('Rad '+str(rads[j])+ ', flx '+str(fl) +' ' + str(fl_blr))
    
            Hal_plots.savefig()
        else:
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
    
    
    np.savetxt(PATH+'Four_Quasars/Results_storage/Growth/'+ID+'Growth_OIII.txt', flux_save)
    
    plt.xlabel('Radius (arcseconds)')
    plt.ylabel('Normalized flux')
    
    #plt.close('all')

plt.show()