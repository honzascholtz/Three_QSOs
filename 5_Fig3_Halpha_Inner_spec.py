#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 18 13:28:28 2018

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


import Fitting_tools as emfit
import IFU_tools_QSO as IFU
import Plotting_tools as emplot

from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.nddata import Cutout2D



#itere = np.array([0])

from matplotlib.colors import LogNorm

Sample = Table.read(PATH+'Four_Quasars/Four_Quasars.fits')
New = Sample

itere= np.array([2,1,0])

ii = 0
jj = 0


f, axes = plt.subplots(3,1, figsize=(5,12))

#itere= np.array([1,4])
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
    ################
    # Spectrum - Halpha
    ################
    
    #ax1 = axes[ii,jj]
    ax1 = axes[jj]
    #plt.setp(ax1.get_yticklabels()[:], visible=False)  
    
    
    if ID=='XID_2028':
        Data = np.loadtxt(PATH+'KMOS_SIN/Results_storage/Spectrums/lid_1565_H_tot.txt')
    else:
        Data = np.loadtxt(PATH+'KMOS_SIN/Results_storage/Spectrums/'+ID+'_H_tot.txt')
        
        Data = np.loadtxt(PATH+'Four_Quasars/Results_storage/Spectrums/'+ID+'_H_inn.txt')
        
    
      
    
    wave = Data[0,:]
    flux = np.ma.array(data=Data[1,:], mask=Data[3,:])
    error = Data[2,:]    
    out ,chi2 = emfit.fitting_Halpha_mul_bkp(wave,flux,error,z, broad=1)
    
    print out.params['Han_fwhm'].value/out.params['Han_center'].value*3e5
    
    z =  (out.params['Han_center'].value*1e4/6562.)-1
    
    
    print out.params['Han_fwhm'].value/out.params['Han_center'].value*3e5
    fl = IFU.flux_measure_ind(out,wave, 'H', use='broad')
       
    print ID, ' flux broad ',  fl
    
    emplot.plotting_Halpha(wave, flux, ax1, out, 'mul',z, title=0)
    
    ax1.set_ylim(-0.1*max(out.eval(x=wave)), 1.2*max(out.eval(x=wave)))
    
    
    ax1.fill_between((6710, 6736), y1 =-0.1*max(out.eval(x=wave)), y2 = 1.2*max(out.eval(x=wave)) , color='grey', alpha=0.4)
    
    ax1.text(6737, max(out.eval(x=wave)), '[SII]', color='firebrick')
    
    ax1.set_xlim(6450, 6750)    
    ax1.set_xlim(6200,7000)
    
    ax1.text(6250, max(out.eval(x=wave)), ID)
    
    
    #if ii==1:
    #    ax1.set_xlabel('Rest wavelength (ang)', fontsize=12)
        
    
    #if jj==0:
    #    ax1.set_ylabel(r'x$10^{-17}$ (erg/s/cm$^2$/ang)', fontsize=12)
    ax1.set_ylabel(r'x$10^{-17}$ (erg/s/cm$^2$/ang)', fontsize=12)
    
    if jj==2:
        ax1.set_xlabel('Rest wavelength (ang)', fontsize=12)
    
    jj+=1   
    #if jj==2:
    #    jj=0
    #    ii=1
    
   
    

plt.savefig(PATH+'Four_Quasars/Graphs/Paper_plots/Fig3_Halpha_spectra.pdf', bbox_inches = 'tight')

plt.show()