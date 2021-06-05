#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 17 10:44:52 2020

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

def create_circular_mask(h, w, center=None, radius=None):

    if center is None: # use the middle of the image
        center = [int(w/2), int(h/2)]
    if radius is None: # use the smallest distance between the center and image walls
        radius = min(center[0], center[1], w-center[0], h-center[1])

    Y, X = np.ogrid[:h, :w]
    dist_from_center = np.sqrt((X - center[0])**2 + (Y-center[1])**2)

    mask = dist_from_center <= radius
    return mask

def create_annular_mask(h, w, center=None, radius=None, radius2=None):

    if center is None: # use the middle of the image
        center = [int(w/2), int(h/2)]
    if radius is None: # use the smallest distance between the center and image walls
        radius = min(center[0], center[1], w-center[0], h-center[1])

    Y, X = np.ogrid[:h, :w]
    dist_from_center = np.sqrt((X - center[0])**2 + (Y-center[1])**2)

    mask = (dist_from_center >= radius) & (dist_from_center <= radius2)
    return mask

def smooth(image,sm):
    
    from astropy.convolution import Gaussian2DKernel
    
    from astropy.convolution import convolve
    
    gauss_kernel = Gaussian2DKernel(sm)

    con_im = convolve(image, gauss_kernel)
    
    #con_im = con_im#*image/image
    
    return con_im  



import Tools_IFU as IFU
import Tools_plotting as emplot
import Tools_fitting as emfit
import Tools_path as ph


from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.nddata import Cutout2D




limx={'2QZJ': [2.15, 2.35, 2.2361, 616]} #
limx['LBQS'] = [2.17,2.24, 2.203, 250]
limx['HB89'] = [2.2,2.3, 2.2560, 500]

lims={'HB89': np.array([-9,14,-11,12])}
lims['2QZJ'] = np.array([-10,10,-10,10])
lims['LBQS'] = np.array([-10,10,-8,12])



Sample = Table.read(PATH+'Four_Quasars/Four_Quasars.fits')
New = Sample

itere= np.array([2,1,0])


f = plt.figure(figsize=(12,10))

column = 0


def forward(x):
    
    return x/(1+z)*1e4

def back(x):
    
    return x*(1+z)/1e4


for i in itere:
    ID = Sample['ID'][i].strip()
    
    if ID=='HB89':
        IDp='HB8903'
    if ID=='2QZJ':
        IDp='2QZJ00'       
    if ID=='LBQS':
        IDp= 'LBQS01'
        
    print( ID) 
    z = Sample['z'][i]
    
    storage = emplot.load_obj('Results_storage/Props/'+ID+'_H')
        
    storage['Cat_entry'] = New[i]
    storage['X-ray ID'] = ID
    
    center = storage['Cube_cent'][1:3]
    
    Data = pyfits.getdata(ph.MyPATH+'Results_storage/Halpha/'+ID+'_Nearest_spaxel_fit_res.fits')
    header = pyfits.getheader(ph.MyPATH+'Results_storage/Halpha/'+ID+'_Nearest_spaxel_fit_res.fits')
    
    n_spixels = header['NAXIS3']
    #dim = [n_ypixels, n_xpixels, n_spixels]
    
    wave = header['CRVAL3'] + (np.arange(n_spixels) - (header['CRPIX3'] - 1.0))*header['CDELT3']
    
    wv_rest = wave/(1+z)*1e4
    
    
    #Data[:, 16,19] = 0
    #Data[:, 14:19,18] = 0
    
    mks = np.loadtxt(ph.MyPATH+ID+'_Skyline_mask.txt')
    rng = limx[ID][3]/3e5*limx[ID][2]
    
# =============================================================================
#    Ring Spectrum
# =============================================================================
    ax3 = f.add_axes( (0.07+column*0.32, 0.07, 0.28, 0.2) )
    center= [center[0], center[1]]
    # Inner region:
    
    h, w = Data.shape[1:3]
    mask = create_annular_mask(h, w, center= center, radius= 3, radius2=8) # Creating the mask 
    mask = np.invert(mask)
    
    mask_in = np.zeros_like(Data)
    
    # Masking the whole cube 
    for j in range(Data.shape[0]):     
        mask_in[j,:,:] = mask
    
    
    Spec_in = np.ma.array(data=Data, mask=mask_in)
    Spec_in = Spec_in.sum(axis=(1,2))
    Spec_in = np.ma.array(data=Spec_in, mask=mks)
    
    fit = np.where( (wave>limx[ID][0]) & (wave<limx[ID][1]) )[0]
    
    wv_plt = wave[fit]
    yplot = Spec_in[fit]
    
    # Converting to velocities
    vels_plt = (wv_plt-limx[ID][2])/limx[ID][2]*3e5
    
    ax3.plot(vels_plt, yplot.data, drawstyle='steps-mid', color='grey', alpha=0.2)
    
    yplotms = yplot.data[np.invert(yplot.mask)]
    wv_pltms = wv_plt[np.invert(yplot.mask)]
    vels_pltms = (wv_pltms-limx[ID][2])/limx[ID][2]*3e5
    
    ax3.plot(vels_pltms, yplotms, drawstyle='steps-mid')
   
    ax3.set_xlim(limx[ID][0],limx[ID][1])
    ax3.set_xlim(-900,900)
       
    ax3.set_ylim(-1,2)
    
    ax3.vlines(0, -1000, 1000, color='k', linestyle='dashed', alpha=0.99)
    ax3.hlines(0, -1000, 1000, color='r', linestyle='dashed', alpha=0.99)
    
    ax3.set_xlabel(r'Velocity shift (km s$^{-1}$)')
    
    if i==0:
        ax3.set_ylim(-0.5,0.5)
        
    if i==2:
        ax3.set_ylabel(r'x$10^{-17}$ (erg/s/cm$^2$/ang)')
    
    
    if i==0:
        ax3.set_ylim(-0.3,0.5)
        ax3.text(-800,0.35, IDp)
        
    else:
        ax3.set_ylim(-1,2)
        
       
        ax3.text(-800,1.5, IDp)
   
    column+=1

     
    ax3.tick_params( direction='in')
    
    
    
    
f.savefig(ph.MyPATH+'Graphs/Paper_plots/Fig7_Residual_spectra.pdf', bbox_inches='tight') 
    


plt.show()    