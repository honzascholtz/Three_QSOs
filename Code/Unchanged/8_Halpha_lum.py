#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 28 13:38:37 2020

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
import scipy.integrate as scpi
nan= float('nan')


import IFU_tools_QSO as IFU
import Fitting_tools as emfit
import Plotting_tools as emplot

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


Halpha_flux_tot = np.zeros(3)
Halpha_flux_nuc = np.zeros(3)
Hbeta_flux_nuc = np.zeros(3)

Av = np.zeros(3)
Halpha_flux_tot_cor = np.zeros(3)

def gauss(x,k,mu,sig):
    expo= -((x-mu)**2)/(sig*sig)
    
    y= k* e**expo
    
    return y

for i in itere:
    ID = Sample['ID'][i].strip()
    z= Sample['z'][i]
# =============================================================================
#     Total Halpha flux
# =============================================================================
    Res = np.loadtxt(PATH+'Four_Quasars/Results_storage/Growth/'+ID+'Growth_Hal.txt')
    
    Flux = Res[1,:]
    
    Halpha_flux_tot[i] = max(Flux)
    
    
# =============================================================================
# Nuclear OIII flux
# =============================================================================
      
    Halpha_flux_nuc[i] = np.loadtxt(PATH+'Four_Quasars/Results_storage/'+ID+'_Halpha_nuclear_flux.txt')
    
    
# =============================================================================
#     Nuclear Hbeta - there are a lot of upper limits
# =============================================================================   
    Save_oiii = np.loadtxt(PATH+'Four_Quasars/Results_storage/Spectrums/'+ID+'_OIII_inn.txt')
    
    wave = Save_oiii[0,:]
    error = Save_oiii[2,:]
    Spectra= np.ma.array(data=Save_oiii[1,:], mask=Save_oiii[3,:])
    
    
    hbw=4861
    if (ID=='LBQS') | (ID=='HB89'):
        hbw=4870.
    out = emfit.fitting_OIII_Hbeta_qso_mul(wave,Spectra, error,z, hbw=hbw)
    
    g, (ax,ax2) = plt.subplots(2,1)
    
    emplot.plotting_OIII_Hbeta(wave, Spectra, ax, out, 'mul',z, title=0)
    
    error_spec =  np.random.normal(0, error[0], len(wave))
    error_spec = np.ma.array(data=error_spec, mask=Save_oiii[3,:])
    
    
    ax2.plot(wave/(1+z)*1e4, error_spec, drawstyle='steps-mid')
    
    ax2.set_xlim(4700, 5100)
    ax2.set_ylim(-0.2, 0.5)
    
    peak_test = np.linspace(0.5*error[0], 5*error[0], 50)  
    
    cent = 4861*(1+z)/1e4
    
    '''
    if (i==0) | (i==1):
        Hbeta_flux_calc = np.array([])
        
        for kl in range(300):      
            for j in range(50):
                error_spec =  np.random.normal(0, error[0], len(wave))
                error_spec = np.ma.array(data=error_spec, mask=Save_oiii[3,:])
                
                flux_test = error_spec + gauss(wave, peak_test[j], cent, 600./3e5*cent)
                
                New = emfit.fitting_Hbeta_sig(wave, flux_test, error,z)                                  
                SNR,chi2 = IFU.SNR_calc(flux_test, wave, New, 'Hbs',z, mul=1)
                                       
                if SNR>3.:
                    break
                            
            Hbeta_flux_calc = np.append(Hbeta_flux_calc, scpi.simps(New.eval_components(x=wave)['Hb_'] , wave)* 1e-13)
    
    
        plt.figure()
        plt.hist(Hbeta_flux_calc, bins=np.linspace(1e-17,2e-16, 20 ))
        
        Hbeta_flux_nuc[i] = np.median(Hbeta_flux_calc)
    
    elif i==2: 
        Hbeta_flux_nuc[i] = float(np.loadtxt(PATH+'Four_Quasars/Results_storage/'+ID+'_Hbeta_nuclear_flux.txt'))
     '''  
        
# =============================================================================
#     Doing actual calculations 
# =============================================================================
    #Av[i] = IFU.Av_calc(Halpha_flux_nuc[i], Hbeta_flux_nuc[i])
    
    #Halpha_flux_tot_cor[i] = IFU.Flux_cor(Halpha_flux_tot[i], Av[i])
    
    


print (Hbeta_flux_nuc)


plt.show()



    



    
    
    
    