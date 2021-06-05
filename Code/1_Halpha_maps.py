#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 16 14:51:35 2017

@author: jscholtz
"""

#importing modules
import numpy as np
import matplotlib.pyplot as plt

from astropy.io import fits as pyfits
from astropy import wcs
from astropy.table import Table, join, vstack, Column
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

global_off= 0.2


def luminosity(z,flux,gamma):
    from astropy.cosmology import FlatLambdaCDM
    import astropy.units as u
    cosmo = FlatLambdaCDM(H0=72 * u.km / u.s / u.Mpc, Om0=0.3)
    d_lum = cosmo.luminosity_distance(z)*(3.08*10**24)
    L= 4*pi*(d_lum**2)*flux*(1+z)**(gamma-2)
    
    return L 

fsz = gst.graph_format()

binning = 'Nearest' 

import Tools_IFU as IFU
import Tools_plotting as emplot
import Tools_fitting as emfit
import Tools_path as ph


plot_it = 0
#0,4

Sample = Table.read(ph.MyPATH+'Four_Quasars.fits')


Hal_size={'HB89': np.array(['0.1', '1'], dtype=str)}
Hal_size['LBQS'] = np.array(['0.25', '1'], dtype=str)
Hal_size['2QZJ'] = np.array(['0.2', '1'], dtype=str)
Hal_size['XID_2028'] = np.array(['0.2', '1'], dtype=str)
   
itere=np.array([0,1,2])

Hal = True
Spat = True

for i in itere:    
    #ID = New['XID'][i]
  

    if Hal== True:
        
        ID = Sample['ID'][i].strip()
        z= Sample['z'][i]
        Hal_band='K_sin'
        
        H_file = Sample['Halpha cube'][i].strip() 
        
        #Loading the Cube
        storage_H = IFU.load_cube(ph.MyPATH+'SINFONI/'+H_file+'.fits', z, ID, 'Sinfoni')
        # Adding additional info to the python dictionary from the QSO catalogue
        storage_H = IFU.add_res(storage_H, Sample[i])
    
        # Masking Emission Region - for collapsing continuum
        storage_H = IFU.mask_emission(storage_H, z)
    
        # Very quick Sky masking - just to collpase the continuum
        storage_H = IFU.mask_sky(storage_H, 1.5)
    
        # Unmasking emission region
        storage_H = IFU.unmask_em(storage_H,z)

    
        # Collapsing the cube to image
        storage_H = IFU.collapse_white(storage_H, plot_it)
        
        # Finding center of the continuum emission 
        storage_H = IFU.find_center(storage_H,'Median_stack_white', plot_it)
    

        size = Hal_size[ID]
        
        rds = float(size[0])
        fl = size[1]
        print ('Radius and extraction region ', rds, fl)
        
        # Extracting spectrum from an aperture
        storage_H = IFU.choose_pixels(storage_H, plot_it, rad= rds , flg = fl)
        
        # Applying astrometry corrections based on GAIA
        storage_H = IFU.astrometry_correction_GAIA(storage_H)
    
        # Stacking sky emission to create proper skyline mask
        storage_H = IFU.stack_sky(storage_H, Hal_band, plot_it, expand=0)
        
        # Collpasing the object spectrum to a single spectrum 
        storage_H = IFU.D1_spectra_collapse(storage_H,z, 'H', plot_it, addsave='_tot')
        
        storage_H = IFU.fitting_collapse_Halpha(storage_H,z, 1, broad=1)
    
        
        Sub_stor= 1
        import Tools_plotting as emplot
        emplot.Summary_plot(storage_H, 'H',z, ID, Sub_stor, prj = '4QSO')
        
        
          
        
        if Spat== True:
            
            storage_H = IFU.Spaxel_fit_sig(storage_H, 'H',1, binning ,broad=1, localised=0 )
    
        
        
        storage_new = IFU.Sub_QSO(storage_H)    
        storage_new = IFU.D1_spectra_collapse(storage_new,z, 'H', 0)
        storage_new= IFU.fitting_collapse_Halpha(storage_new,z, 1, broad=0)
        
        f = plt.figure(figsize=(8,4))
        
        ax1 = f.add_subplot('121')
        
        ax1.tick_params(direction='in')       
        wave = storage_new['obs_wave'].copy()
        flux = storage_new['1D_spectrum'].copy()    
        out = storage_new['1D_fit_Halpha_mul'] 
        emplot.plotting_Halpha(wave, flux, ax1, out,'mul',z)
        
        mns = np.ma.min(flux).copy()
        
        
        
        wave = storage_H['obs_wave'].copy()/(1+z)*1e4
        flux = storage_H['1D_spectrum'].copy()    
        
        ax1.plot(wave,flux, drawstyle='steps-mid')
        
        try:
            ax1.set_ylim(1.*mns, 1.1*np.ma.max(flux))
            
        except:
            ax1.set_ylim( -0.1, 0.4)
        
        
        header_cube = storage_H['header']
        cube_wcs= wcs.WCS(header_cube).celestial 
        
        ax = f.add_subplot('122', projection=cube_wcs)
 
        out = storage_H['1D_fit_Halpha_mul']
        wv_obs = storage_H['obs_wave'].copy()
        
        cent = out.params['Han_center'].value
        sig = out.params['Han_sigma'].value
        
        Cube_ready = storage_new['flux'].data
        Cube_ready = np.ma.array(data=Cube_ready, mask=storage_new['sky_clipped'])*1e-13
        
        use = np.where( (wv_obs<cent+(1*sig)) & (wv_obs>cent-(1*sig)))[0]
        
        muls = len(use)*(wv_obs[1]-wv_obs[0])
        
        print (muls)
        
        hal = np.ma.sum(Cube_ready[use,:,:], axis=(0)) *muls
        
        
        if ID== '2QZJ':
            ax.imshow(hal, origin='low', vmax=8e-18, vmin=0)
        
        else:
            ax.imshow(hal, origin='low')
        
        
        
        emplot.new_ALMA_contour_plot(storage_H, ax, prj='4QSO') 
        
        f.savefig(ph.MyPATH+'Graphs/'+ID+'_Sub_QSO_spec_map.pdf')
        
        
        np.savetxt(ph.MyPATH+'Results_storage/Sub_qso_map_hal/'+ID+'.txt', hal)
        
        
        
    
plt.show()
    

        
            




