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

import Graph_setup as gst 
import Tools_IFU as IFU
import Tools_plotting as emplot
import Tools_fitting as emfit
import Tools_path as ph


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




from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.nddata import Cutout2D




limx={'2QZJ': [2.15, 2.35, 2.2361, 616]} #
limx['LBQS'] = [2.17,2.24, 2.203, 250]
limx['HB89'] = [2.2,2.3, 2.2560, 500]

lims={'HB89': np.array([-9,14,-11,12])}
lims['2QZJ'] = np.array([-10,10,-10,10])
lims['LBQS'] = np.array([-10,10,-8,12])



Sample = Table.read(ph.MyPATH+'Four_Quasars.fits')
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
    
    Data = pyfits.getdata(ph.MyPATH+'Results_storage/Hal_fitting/'+ID+'_Nearest_spaxel_fit_res.fits')
    header = pyfits.getheader(ph.MyPATH+'Results_storage/Hal_fitting/'+ID+'_Nearest_spaxel_fit_res.fits')
    
    n_spixels = header['NAXIS3']
    #dim = [n_ypixels, n_xpixels, n_spixels]
    
    wave = header['CRVAL3'] + (np.arange(n_spixels) - (header['CRPIX3'] - 1.0))*header['CDELT3']
    
    wv_rest = wave/(1+z)*1e4
    
    
    #Data[:, 16,19] = 0
    #Data[:, 14:19,18] = 0
    
    mks = np.loadtxt(ph.MyPATH+ID+'_Skyline_mask.txt')
# =============================================================================
#     Plotting setup
# =============================================================================
    IFU_wcs= wcs.WCS(header).celestial
    
    ax1 = f.add_axes( (0.08+column*0.32, 0.69, 0.25, 0.3), projection=IFU_wcs )
    
    cmap = plt.get_cmap('cividis')
    cmap.set_bad(color='black')
    
# =============================================================================
#   Imaging the location 
# =============================================================================
    #collapse = np.where((wv_rest<6562+5) & (wv_rest>6562-5))[0]
    
    rng = limx[ID][3]/3e5*limx[ID][2]
    
    collapse = np.where((wave<limx[ID][2]+rng/2) & (wave>limx[ID][2]-rng/2))[0]
     
    Res_image = Data[collapse,:,:].sum(axis=0)
    
    #Res_image = smooth(Res_image, 1)
    
    Res_image[np.where(Res_image==0)] = np.nan
      
    rms = np.nanstd(Res_image- np.nanmean(Res_image))
    
    
    ax1.imshow(Res_image,vmin=-2*rms, vmax=3*rms,cmap=cmap,   origin='low')
    
    ax1.contour(Res_image, levels=(2*rms, 3*rms, 5*rms), linestyles='dashed', colors='red')
    
    lim = lims[ID]
    coor= center
    ax1.set_xlim(coor[0]+lim[0],coor[0]+lim[1])
    ax1.set_ylim(coor[1]+lim[2],coor[1]+lim[3])
    
    deg_per_pix = header['CDELT2']
    arc_per_pix = deg_per_pix*3600
    
    dx = lim[1] - lim[0]
    dy = lim[3] - lim[2]
    x = coor[0]+lim[0]
    y = coor[1]+lim[2]
        
    
    lim_sc = lim*arc_per_pix
    
    emplot.overide_axes_labels(f,ax1,(lim_sc), labelx=1, labely=1,tickin=1)
    
    ax1.text(dx*0.05+x, dy*0.9+y, IDp, color='white')
# =============================================================================
#   Nuclear spectrum
# =============================================================================
    ax2 = f.add_axes( (0.06+column*0.32, 0.4, 0.28, 0.2) )
    center= [center[0], center[1]]
    # Inner region:
    rad = 3
    
    h, w = Data.shape[1:3]
    mask = create_circular_mask(h, w, center= center, radius= rad)
    mask = np.invert(mask)
    
    mask_in = np.zeros_like(Data)
    
    for j in range(Data.shape[0]):
        
        mask_in[j,:,:] = mask
 
    
    Spec_in = np.ma.array(data=Data, mask=mask_in)
    Spec_in = Spec_in.sum(axis=(1,2))
    Spec_in = np.ma.array(data=Spec_in, mask=mks)
    
    
    fit = np.where( (wave>limx[ID][0]) & (wave<limx[ID][1]) )[0]
    
    
    fitx = wave[fit]
    fity = Spec_in[fit]    
    
    p_coeff = np.polyfit( fitx, fity, 2 ) 
    p = np.poly1d( p_coeff ) 
    
    
    
    wv_plt = wave[fit]
    yplot = Spec_in[fit]#-p(fitx)
    

    ax2.plot(wv_plt, yplot.data, drawstyle='steps-mid', color='grey', alpha=0.2)
    
    yplotms = yplot.data[np.invert(yplot.mask)]
    wv_pltms = wv_plt[np.invert(yplot.mask)]
    
    ax2.plot(wv_pltms, yplotms, drawstyle='steps-mid')
    #ax2.plot(fitx, p(fitx))
    #ax2.plot(wave, y)
    
    ax2.vlines(limx[ID][2], -2,2)
    
    ax2.vlines(limx[ID][2]-rng/2, -2,2, linestyle='dashed', color='red')
    ax2.vlines(limx[ID][2]+rng/2, -2,2, linestyle='dashed', color='red')
    
    
    ax2.set_xlim(limx[ID][0],limx[ID][1])
    #ax2.set_xlim(6400, 6686)
    #print(limx[ID]/(1+z)*1e4)
    
    if i==0:
        ax2.set_ylim(-0.3,0.5)
        
        ax2.text(limx[ID][0]+0.01,0.35, 'Nuclear Spectrum')
        ax2.text(limx[ID][0]+0.005,0.25, 'r<0.5')
        
    else:
        ax2.set_ylim(-1,2)
        ax2.text(limx[ID][0]+0.005,1.5, 'Nuclear Spectrum')
        
        ax2.text(limx[ID][0]+0.005,1.1, 'r<0.5')
        

        
    ax2.set_xlabel(r'Observed Wavelength ($\mu$m)')    
    axt = ax2.secondary_xaxis('top', functions=(forward, back))
    axt.set_xlabel('Rest wavelength (A)')
    axt.tick_params( direction='in')
    #axt.set_xlim(limx[ID][0]/(1+z)*1e4,limx[ID][1]/(1+z)*1e4)  
# =============================================================================
#    Ring Spectrum
# =============================================================================
    ax3 = f.add_axes( (0.06+column*0.32, 0.07, 0.28, 0.2) )
    center= [center[0], center[1]]
    # Inner region:
    rad = 3
    
    h, w = Data.shape[1:3]
    mask = create_annular_mask(h, w, center= center, radius= rad, radius2=8)
    mask = np.invert(mask)
    
    mask_in = np.zeros_like(Data)
    
    for j in range(Data.shape[0]):
        
        mask_in[j,:,:] = mask
    
    
    
    Spec_in = np.ma.array(data=Data, mask=mask_in)
    Spec_in = Spec_in.sum(axis=(1,2))
    Spec_in = np.ma.array(data=Spec_in, mask=mks)
    
    fit = np.where( (wave>limx[ID][0]) & (wave<limx[ID][1]) )[0]
      
    
    fitx = wave[fit]
    fity = Spec_in[fit]    
    
    p_coeff = np.polyfit( fitx, fity, 2 ) 
    p = np.poly1d( p_coeff ) 
    
    wv_plt = wave[fit]
    yplot = Spec_in[fit]-p(fitx)
    
    ax3.plot(wv_plt, yplot.data, drawstyle='steps-mid', color='grey', alpha=0.2)
    
    yplotms = yplot.data[np.invert(yplot.mask)]
    wv_pltms = wv_plt[np.invert(yplot.mask)]
    
    ax3.plot(wv_pltms, yplotms, drawstyle='steps-mid')
    
    #ax3.vlines(limx[ID][2], -2,2)
    
    ax3.vlines(limx[ID][2]-rng/2, -2,2, linestyle='dashed', color='red')
    ax3.vlines(limx[ID][2]+rng/2, -2,2, linestyle='dashed', color='red')
    
    
    ax3.set_xlim(limx[ID][0],limx[ID][1])
       
    ax3.set_ylim(-1,2)
    
    ax3.set_xlabel(r'Observed Wavelength ($\mu$m)')
    
    if i==0:
        ax3.set_ylim(-0.5,0.5)
        
        
    if i==2:
        
        ax2.set_ylabel('Flux density')
        ax3.set_ylabel('Flux density')
    
    
    if i==0:
        ax3.set_ylim(-0.3,0.5)
        
        ax3.text(limx[ID][0]+0.01,0.35, 'Ring Spectrum')
        ax3.text(limx[ID][0]+0.005,0.25, '0.5<r<0.8')
    else:
        ax3.set_ylim(-1,2)
        ax3.text(limx[ID][0]+0.005,1.5, 'Ring Spectrum')
        
        ax3.text(limx[ID][0]+0.005,1.1, '0.5<r<0.8')
    
    axt = ax3.secondary_xaxis('top', functions=(forward, back))
    axt.set_xlabel('Rest wavelength (A)')
    axt.tick_params( direction='in')
    column+=1
# =============================================================================
# Testing  
# =============================================================================
    
    
    Image_narrow_old = pyfits.getdata(PATH+'KMOS_SIN/Results_storage/Halpha/'+ID+'_Nearest_spaxel_fit.fits')
    
    Hal_map = Image_narrow_old[5]
    
    mx = np.nanmax( Hal_map)
    
    ax1.contour(Hal_map, levels = (0.1*mx,0.68*mx, 0.87*mx, 0.9*mx), colors='white' )
     


    ax2.tick_params( direction='in')
    ax3.tick_params( direction='in')
    
    
    
    
#f.savefig(ph.MyPATH+'Graphs/Paper_plots/Residual_emission.pdf') 
    
''' 
# =============================================================================
# =============================================================================
# =============================================================================
# =============================================================================
# =============================================================================
# =============================================================================
# =============================================================================
# =============================================================================
# =============================================================================
# # # # # # # # #     
# =============================================================================
# =============================================================================
# =============================================================================
# =============================================================================
# =============================================================================
# =============================================================================
# =============================================================================
# =============================================================================
# =============================================================================
itere = np.array([2])
limx={'2QZJ': [2.15, 2.35, 2.2361, 616]} #
limx['LBQS'] = [2.17,2.24, 2.203, 250]
limx['HB89'] = [2.2,2.3, 2.2505, 700]


for i in itere:
    ID = Sample['ID'][i].strip()
    
    print( ID) 
    z = Sample['z'][i]
    
    storage = emplot.load_obj('KMOS_SIN/Results_storage/Props/'+ID+'_H')
        
    storage['Cat_entry'] = New[i]
    storage['X-ray ID'] = ID
    
    center = storage['Cube_cent'][1:3]
    
    Data = pyfits.getdata(PATH+'Four_Quasars/Results_storage/Hal_fitting/'+ID+'_Nearest_spaxel_fit_res_2L.fits')
    header = pyfits.getheader(PATH+'Four_Quasars/Results_storage/Hal_fitting/'+ID+'_Nearest_spaxel_fit_res_2L.fits')
    
    n_spixels = header['NAXIS3']
    #dim = [n_ypixels, n_xpixels, n_spixels]
    
    wave = header['CRVAL3'] + (np.arange(n_spixels) - (header['CRPIX3'] - 1.0))*header['CDELT3']
    
    wv_rest = wave/(1+z)*1e4
    
    
    #Data[:, 16,19] = 0
    #Data[:, 14:19,18] = 0
    
    mks = np.loadtxt(PATH+'Four_quasars/'+ID+'_Skyline_mask.txt')
    
    f, (ax1, ax2,ax3) = plt.subplots(3, figsize=(4,8))
# =============================================================================
#     Plotting setup
# =============================================================================
    
    
    cmap = plt.get_cmap('cividis')
    cmap.set_bad(color='black')
    
# =============================================================================
#   Imaging the location 
# =============================================================================
    #collapse = np.where((wv_rest<6562+5) & (wv_rest>6562-5))[0]
    
    rng = limx[ID][3]/3e5*limx[ID][2]
    
    collapse = np.where((wave<limx[ID][2]+rng/2) & (wave>limx[ID][2]-rng/2))[0]
     
    Res_image = Data[collapse,:,:].sum(axis=0)
    
    #Res_image = smooth(Res_image, 1)
    
    Res_image[np.where(Res_image==0)] = np.nan
      
    rms = np.nanstd(Res_image- np.nanmean(Res_image))
    
    
    ax1.imshow(Res_image,vmin=-2*rms, vmax=3*rms,cmap=cmap,   origin='low')
    
    ax1.contour(Res_image, levels=(2*rms, 3*rms, 5*rms), linestyles='dashed', colors='red')
    
       
    lim = lims[ID]
    coor= center
    ax1.set_xlim(coor[0]+lim[0],coor[0]+lim[1])
    ax1.set_ylim(coor[1]+lim[2],coor[1]+lim[3])
    
    deg_per_pix = header['CDELT2']
    arc_per_pix = deg_per_pix*3600
      
    
# =============================================================================
#   Nuclear spectrum
# =============================================================================
    
    center= [center[0], center[1]]
    # Inner region:
    rad = 3
    
    h, w = Data.shape[1:3]
    mask = create_circular_mask(h, w, center= center, radius= rad)
    mask = np.invert(mask)
    
    mask_in = np.zeros_like(Data)
    
    for j in range(Data.shape[0]):
        
        mask_in[j,:,:] = mask
 
    
    Spec_in = np.ma.array(data=Data, mask=mask_in)
    Spec_in = Spec_in.sum(axis=(1,2))
    Spec_in = np.ma.array(data=Spec_in, mask=mks)
    
    
    fit = np.where( (wave>limx[ID][0]) & (wave<limx[ID][1]) )[0]
    
    
    fitx = wave[fit]
    fity = Spec_in[fit]    
    
    p_coeff = np.polyfit( fitx, fity, 2 ) 
    p = np.poly1d( p_coeff ) 
    
    
    
    wv_plt = wave[fit]
    yplot = Spec_in[fit]-p(fitx)
    

    ax2.plot(wv_plt, yplot.data, drawstyle='steps-mid', color='grey', alpha=0.2)
    
    yplotms = yplot.data[np.invert(yplot.mask)]
    wv_pltms = wv_plt[np.invert(yplot.mask)]
    
    ax2.plot(wv_pltms, yplotms, drawstyle='steps-mid')
    #ax2.plot(fitx, p(fitx))
    #ax2.plot(wave, y)
    
    ax2.vlines(limx[ID][2], -1,2)
    
    ax2.vlines(limx[ID][2]-rng/2, -2,2, linestyle='dashed', color='red')
    ax2.vlines(limx[ID][2]+rng/2, -2,2, linestyle='dashed', color='red')
    
    
    ax2.set_xlim(limx[ID][0],limx[ID][1])
    
    ax2.set_ylim(-1,2)
    
           
# =============================================================================
#    Ring Spectrum
# =============================================================================
    
    center= [center[0], center[1]]
    # Inner region:
    rad = 3
    
    h, w = Data.shape[1:3]
    mask = create_annular_mask(h, w, center= center, radius= rad, radius2=8)
    mask = np.invert(mask)
    
    mask_in = np.zeros_like(Data)
    
    for j in range(Data.shape[0]):
        
        mask_in[j,:,:] = mask
    
    
    
    Spec_in = np.ma.array(data=Data, mask=mask_in)
    Spec_in = Spec_in.sum(axis=(1,2))
    Spec_in = np.ma.array(data=Spec_in, mask=mks)
    
    fit = np.where( (wave>limx[ID][0]) & (wave<limx[ID][1]) )[0]
      
    
    fitx = wave[fit]
    fity = Spec_in[fit]    
    
    p_coeff = np.polyfit( fitx, fity, 2 ) 
    p = np.poly1d( p_coeff ) 
    
    wv_plt = wave[fit]
    yplot = Spec_in[fit]-p(fitx)
    
    ax3.plot(wv_plt, yplot.data, drawstyle='steps-mid', color='grey', alpha=0.2)
    
    yplotms = yplot.data[np.invert(yplot.mask)]
    wv_pltms = wv_plt[np.invert(yplot.mask)]
    
    ax3.plot(wv_pltms, yplotms, drawstyle='steps-mid')
    
    #ax3.vlines(limx[ID][2], -2,2)
    
    ax3.vlines(limx[ID][2]-rng/2, -2,2, linestyle='dashed', color='red')
    ax3.vlines(limx[ID][2]+rng/2, -2,2, linestyle='dashed', color='red')
    
    
    ax3.set_xlim(limx[ID][0],limx[ID][1])
       
    ax3.set_ylim(-1,2)
    
    
    ax3.set_xlabel(r'Observed Wavelength ($\mu$m)')
    
    
    column+=1
# =============================================================================
# Testing  
# =============================================================================
    
    
    Image_narrow_old = pyfits.getdata(PATH+'KMOS_SIN/Results_storage/Halpha/'+ID+'_Nearest_spaxel_fit.fits')
    
    Hal_map = Image_narrow_old[5]
    
    mx = np.nanmax( Hal_map)
    
    ax1.contour(Hal_map, levels = (0.1*mx,0.68*mx, 0.87*mx, 0.9*mx), colors='white' )
 '''  

plt.show()    