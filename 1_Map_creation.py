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

# Setting my style of plotting in pyplot
import Graph_setup as gst
plt.close('all')
# Implementing my style of plotting and setting default fontsize called fsz
fsz = gst.graph_format()

# Setting few general variables
nan= float('nan')
pi= np.pi
e= np.e


c= 3.*10**8
h= 6.62*10**-34
k= 1.38*10**-23

Ken98= (4.5*10**-44)
Conversion2Chabrier=1.7 # Also Madau
Calzetti12= 2.8*10**-44
arrow = u'$\u2193$' 


# Modifiying this path
PATH='/Users/jansen/Google Drive/Astro/'




# Importing other modules necessary to complete the IFU analyses
import IFU_tools_QSO as IFU
import Plotting_tools as emplot
import Fitting_tools as emfit

# This command bins the spaxel in 0.2 arcsecond radius. If you want to switch of the binning replace with 'Individual'
binning = 'Nearest'

# This switches on all check plots - replace with 0 to switch it off
plot_it = 1
#0,4

# Import catalogue with all the info.
Sample = Table.read(PATH+'Four_Quasars/Four_Quasars.fits')

# Modifying the name for the Halpha cube of 2QZJ
Sample['Halpha cube'][0] = '2QZ0028-28'

Spat = False

# Disctionary full of values for the inner 5 kpc aperture
Hal_size={'HB89': np.array(['0.2', '1'], dtype=str)}
Hal_size['LBQS'] = np.array(['0.2', '1'], dtype=str)
Hal_size['2QZJ'] = np.array(['0.2', '1'], dtype=str)
Hal_size['XID_2028'] = np.array(['0.2', '1'], dtype=str)

Hfdge={'HB89': 21}
Hfdge['LBQS'] = 30
Hfdge['2QZJ'] = 27
   
itere=np.array([2])
for i in itere:  
    
    # Extracts main information from the catalogue
    ID = Sample['ID'][i].strip()
    z= Sample['z'][i]
    Hal_band='K_sin'
    
    H_file = Sample['Halpha cube'][i].strip() 
    
    # Loads the cube and creates the disctionary called storage. Following info will be 
    storage_H = IFU.load_cube(PATH+'Four_Quasars/SINFONI/'+H_file+'.fits', z, ID, 'Sinfoni')

    # Adds soem information from the table to the storage
    storage_H = IFU.add_res(storage_H, Sample[i])

    # Masking Emission Region
    storage_H = IFU.mask_emission(storage_H, z)

    # Very quick Sky masking - Just a quick skyline masking to have better continuum image. More theral job will be done later
    storage_H = IFU.mask_sky(storage_H, 1.5)

    # Unmasking emission region
    storage_H = IFU.unmask_em(storage_H,z)


    # Collapsing the cube to create a continuum image
    storage_H = IFU.collapse_white(storage_H, plot_it)
    
    # Finding the center of the continuum image
    storage_H = IFU.find_center(storage_H,'Median_stack_white', plot_it)


    size = Hal_size[ID]
    
    rds = float(size[0])
    fl = size[1]
    print 'Radius and extraction region ', rds, fl
    # Creating a mask which indicate the aperture
    storage_H = IFU.choose_pixels(storage_H, plot_it, rad= rds , flg = fl)
    
    # Correcting the astrometry based on GAIA DR2 catalogue
    storage_H = IFU.astrometry_correction_GAIA(storage_H)
    

    # Stacks all the spaxel outside the object to create a sky spectrum. We then mask everything outside of 1sigma.
    storage_H = IFU.stack_sky(storage_H, Hal_band, plot_it, expand=0)
    
    # Collapses the cube in the aperture specified in Choose_pixels yp create a object spectrum
    storage_H = IFU.D1_spectra_collapse(storage_H,z, 'H', plot_it, addsave='_tot')
    
    # Fitting the collpased object spectrum with models from Fitting_Tools
    storage_H = IFU.fitting_collapse_Halpha(storage_H,z, plot_it, broad=1)

    # Makking summary plots of the fitting and skyline plotting
    Sub_stor= 1
    emplot.Summary_plot(storage_H, 'H',z, ID, Sub_stor, prj = '4QSO')
    
    
          
    # Fitting the speaxel spectra to create a map
    #storage_H = IFU.Spaxel_fit_sig(storage_H, 'H',1, binning ,broad=1, localised=1 )
    # Displaying all the maps created by the Spaxel_fit_sig  
    #emplot.Spax_fit_plot_sig(storage_H, 'H', z, ID, binning, prj='4QSO')
      
    # Subtracting the QSO
    storage_new = IFU.Sub_QSO(storage_H) 
    # Collapsing the subtracted cube to create a new BLR and continuum subtracted spectrum
    
    storage_new = IFU.D1_spectra_collapse(storage_new,z, 'H', 0)
    # Fitting the new spectrum 
    storage_new = IFU.fitting_collapse_Halpha(storage_new,z, 1, broad=0)
    
# =============================================================================
#     Plotting the new BLR and continuum subtracted map
# =============================================================================
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
    
    print muls
    
    fudge = Hfdge[ID]
    hal = np.ma.sum(Cube_ready[use,:,:], axis=(0)) *muls/fudge
    
    
    if ID== '2QZJ':
        ax.imshow(hal, origin='low', vmax=8e-18, vmin=0)
    
    else:
        ax.imshow(hal, origin='low')
    
    
    
    emplot.new_ALMA_contour_plot(storage_H, ax, prj='4QSO') 
    
    #f.savefig(PATH+'Four_Quasars/Graphs/'+ID+'_Sub_QSO_spec_map.pdf')
    
    
    
    f, ax = plt.subplots(1)
    
    ax.imshow(storage_new['BLR_map'], origin='low')
    
    
    
    
    #np.savetxt(PATH+'Four_Quasars/Results_storage/Sub_qso_map_hal/'+ID+'.txt', hal)
        
    
        
    
plt.show()
    

        
            




