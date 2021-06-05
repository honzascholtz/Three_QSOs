#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 14 14:39:26 2020

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

import Tools_path as ph

fsz = gst.graph_format()

###################################
# Need to recalculate things from FWHM to W80 parameters to plot in Figure 1
###################################
# Loading catalogue
#SDSS = Table.read(ph.MyPATH+'Catalogues/ALPAKA_v1.fits')
SDSS = Table.read(ph.MyPATH+'Catalogues/ALPAKA_XMM.fits')


# Setting up the loot to loop over
ind = range(len(SDSS))

# Creating a new table 
New_tab = Table()

# Arrays to store new quantities in 
OIII_FWHM = np.zeros(len(ind))
OIII_W80 = np.zeros(len(ind))
OIII_LUM = np.zeros(len(ind))

import scipy.integrate as scpi

def gauss(x,k,mu,sig):
    expo= -((x-mu)**2)/(sig*sig)
    
    y= k* e**expo
    
    return y

def find_nearest_ind(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx


for i in ind:
    
    if (i % 100 == 0 ) == True:
        print (i)
    
    dat = SDSS[i]
    
    N = 1000
    # Creating a model spectrum based on the paramters from the table
    cent = cntd =  dat['OIII_5007_WAV'] 
    wvs = x = np.linspace( cntd-(3000./3e5*cntd), cntd+(3000./3e5*cntd) , N)
     
    siga = dat['OIII_5007_FWHM']/2.235/3e5*cntd
    sigb = dat['OIII_5007B_FWHM']/2.235/3e5*cntd
    
    OIII_LUM[i] = (dat['OIII_5007_LUM'] + dat['OIII_5007B_LUM'])  
    y = gauss(x , dat['OIII_5007_NORM'], dat['OIII_5007_WAV'], siga) + gauss(x , dat['OIII_5007B_NORM'], dat['OIII_5007B_WAV'], sigb)
    
    #Creating  an array to store the cumulative flux from the
    Int = np.zeros(N-1)
        
    for j in range(N-1):
        
        Int[j] = scpi.simps(y[:j+1], wvs[:j+1]) * 1e-13
                
    Int = Int/max(Int)

    # Finding the array indicaes corresponding to 10, 90 and 50 % of the flux
    ind10 = find_nearest_ind(Int, 0.1)        
    ind90 = find_nearest_ind(Int, 0.9)
    ind50 = find_nearest_ind(Int, 0.5)     
    
    # Finding the wavelengths correspoding to the flux values
    wv10 = wvs[ind10]
    wv90 = wvs[ind90]
    wv50 = wvs[ind50]
    
    # Converting them to velocities
    v10 = (wv10-cent)/cent*3e5
    v90 = (wv90-cent)/cent*3e5
    v50 = (wv50-cent)/cent*3e5
     
    # Calculating the W80 
    w80 = v90-v10
    
    OIII_W80[i] = w80
    
    OIII_FWHM[i] = (dat['OIII_5007_FWHM']*dat['OIII_5007_LUM'] + dat['OIII_5007B_FWHM']*dat['OIII_5007B_LUM'])/OIII_LUM[i]
    
    
New_tab['RA'] = SDSS['RA_1']
New_tab['DEC'] = SDSS['DEC_1']
New_tab['Z'] = SDSS['Z']

New_tab['OIII_W80'] = OIII_W80
New_tab['OIII_LUM'] = OIII_LUM
New_tab['OIII_FWHM'] = OIII_FWHM


New_tab.write(ph.MyPATH+'Catalogues/ALPAKA_W80_XMM.fits', overwrite=True)

plt.show()
    
    
    
