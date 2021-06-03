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


import Graph_setup as gst 
import Tools_IFU as IFU
import Tools_plotting as emplot
import Tools_fitting as emfit
import Tools_path as ph

from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.nddata import Cutout2D



#itere = np.array([0])

from matplotlib.colors import LogNorm

Sample = Table.read(ph.MyPATH+'Four_Quasars.fits')
New = Sample

itere= np.array([2,1,0])

ii = 0
jj = 0


#f, axes = plt.subplots(3,1, figsize=(5,12))

f = plt.figure(figsize=(8,9.6))

k=0
#itere= np.array([1,4])
for i in itere:
    ID = Sample['ID'][i].strip()
    
    z = Sample['z'][i]
    

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
    ax1 = f.add_axes((0.075, 0.695-jj*0.325 ,0.425 ,0.3))
    #ax1 = f.add_axes((0.0625, 0.745-jj*0.325 ,0.425 ,0.25))
    #axres = f.add_axes((0.0625, 0.695-jj*0.325 ,0.425 ,0.05))
    #plt.setp(ax1.get_yticklabels()[:], visible=False)  
    #axres.tick_params(direction='in')
    
    #plt.setp(ax1.get_xticklabels()[:], visible=False) 
    ax1.tick_params(direction='in')
    if i==0:
        xwaaaa=0
    else:
        plt.setp(ax1.get_xticklabels()[:], visible=False) 
    
    
    Data = np.loadtxt(ph.MyPATH+'Results_storage/Spectrums/'+ID+'_H_inn.txt')
    
    wave = Data[0,:]
    flux = np.ma.array(data=Data[1,:], mask=Data[3,:])
    error = Data[2,:]    

    
    #out ,chi2 = emfit.fitting_Halpha_mul_bkp_only(wave,flux,error,z, broad=1)
    
    #out ,chi2 = emfit.fitting_Halpha_mul_LBQS(wave,flux,error,z, broad=1)
    out, chi2 = emfit.fitting_Halpha_mul_testing(wave,flux,error,z, broad=1)
    
    z =  (out.params['Han_center'].value*1e4/6562.)-1
    
    
    emplot.plotting_Halpha(wave, flux, ax1, out, 'mul',z, title=0)
    
    ax1.set_ylim(-0.1*max(out.eval(x=wave)), 1.2*max(out.eval(x=wave)))
    
    
    ax1.fill_between((6710, 6736), y1 =-0.1*max(out.eval(x=wave)), y2 = 1.2*max(out.eval(x=wave)) , color='grey', alpha=0.4)
    
    ax1.text(6737, max(out.eval(x=wave)), '[SII]', color='firebrick')
    
    ax1.set_xlim(6300,6900)
    
    ax1.text(6320, 1.1*max(out.eval(x=wave)), IDp+r' - H$\alpha$', fontsize=15)
    
    ax1.plot(wave/(1+z)*1e4, out.eval_components(x=wave)['X_'], 'k--' )
    
   
    ax1.set_ylabel(r'x$10^{-17}$ (erg/s/cm$^2$/ang)', fontsize=15)
    
    if jj==2:
        ax1.set_xlabel('Rest wavelength (ang)', fontsize=15)
        
    
    wv_rest = wave/(1+z)*1e4
      
    flux_pl = flux.data[np.invert(flux.mask)]
    wv_rst_sc= wv_rest[np.invert(flux.mask)]
    
    fit_loc = np.where((wv_rst_sc>6200)&(wv_rst_sc<7000))[0]
    
    ################
    # Spectrum - OIII
    ################
   
    if ID=='2QZJ':
        z= 2.4064
         
    ax1 = f.add_axes((0.55, 0.695-jj*0.325 ,0.425 ,0.3))
    #ax1 = f.add_axes((0.55, 0.745-jj*0.325 ,0.425 ,0.25))
    #axres = f.add_axes((0.55, 0.695-jj*0.325 ,0.425 ,0.05))
    #plt.setp(ax1.get_yticklabels()[:], visible=False)  
    
    ax1.tick_params(direction='in')
    if i==0:
        xwaaaa=0
    else:
        plt.setp(ax1.get_xticklabels()[:], visible=False) 
    
    Saves = emplot.load_obj('Results_storage/Spectrums/Regions/'+ID+'_OIII_'+str(1)+'_'+str(1))
    
    wave = Saves['wave']/(1+z)*1e4
    Spec = Saves['data']
    
    fit = Saves['total']
    Hbw = Saves['Hbw']
    Hbn = Saves['Hbn']
    
    o3rn = Saves['o3rn']
    o3rw = Saves['o3rw']
    
    o3bn = Saves['o3bn']
    o3bw = Saves['o3bw']
    
    
    ax1.plot(wave, Spec.data,  drawstyle='steps-mid', linestyle='solid', color='grey', alpha=0.2)
    
    flux = Spec.data[np.invert(Spec.mask)]
    wv_rst_sc= wave[np.invert(Spec.mask)]
    
    ax1.plot(wv_rst_sc, flux, drawstyle='steps-mid', linestyle='solid')
    
    ax1.plot(wave, (o3bn), 'green', linestyle='dashed')
    ax1.plot(wave, (o3rn), 'green', linestyle='dashed')
    
    ax1.plot(wave, (o3rw), 'blue', linestyle='dashed')
    ax1.plot(wave, (o3bw), 'blue', linestyle='dashed')
    
    
    ax1.plot(wave, Hbn, 'orange', linestyle='dashdot')
    ax1.plot(wave, Hbw, 'orange', linestyle='dotted')
    
    ax1.plot(wave, fit, 'red', linestyle='dashed')
    
    ax1.set_ylim(-0.1*max(fit), 1.1*max(fit))
    
    ax1.set_xlim(4700,5100)
    
    ax1.text(4710, max(fit), IDp+r' - [OIII] + H$\beta$', fontsize=15)
    
    if jj==2:
        ax1.set_xlabel('Rest wavelength (ang)', fontsize=15)
    
    jj+=1   
    
   
    

plt.savefig(ph.MyPATH+'Graphs/Paper_plots/Fig3_Halpha_OIII_spectra.pdf')#, bbox_inches = 'tight')

plt.show()