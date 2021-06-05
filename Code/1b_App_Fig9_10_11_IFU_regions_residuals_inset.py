#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 16 14:51:35 2017

@author: jscholtz
"""

#importing modules
import numpy as np
import matplotlib.pyplot as plt ; plt.ioff()

from astropy.io import fits as pyfits
from astropy import wcs
from astropy.table import Table, join, vstack, Column


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

PATH='/Users/jansen/Google Drive/Astro/'


def smooth(image,sm):
    
    from astropy.convolution import Gaussian2DKernel
    
    from astropy.convolution import convolve
    
    gauss_kernel = Gaussian2DKernel(sm)

    con_im = convolve(image, gauss_kernel)
    
    #con_im = con_im#*image/image
    
    return con_im  
 
binning = 'Nearest'

import Graph_setup as gst 
import Tools_IFU as IFU
import Tools_plotting as emplot
import Tools_fitting as emfit
import Tools_path as ph

fsz = gst.graph_format()

plot_it = 0
#0,4

Sample = Table.read(ph.MyPATH+'Four_Quasars.fits')
Sample['Halpha cube'][0] = '2QZ0028-28'


OIII_size={'HB89': np.array(['1.2', '1'], dtype=str)}
OIII_size['LBQS'] = np.array(['1.2', '1'], dtype=str)
OIII_size['2QZJ'] = np.array(['1.2', '1'], dtype=str)

Hal_size={'HB89': np.array(['0.2', '1'], dtype=str)}
Hal_size['LBQS'] = np.array(['0.2', '1'], dtype=str)
Hal_size['2QZJ'] = np.array(['0.2', '1'], dtype=str)

lims={'HB89': np.array([-9,14,-11,12])}
lims['2QZJ'] = np.array([-15,15,-15,15])
lims['LBQS'] = np.array([-10,10,-8,12])

IFU_psfs = {'HB89': np.array([-1.1,-0.6,-0.08])/3600}
IFU_psfs['2QZJ'] = np.array([-0.8,-0.7, -0.1 ])/3600
IFU_psfs['LBQS'] = np.array([ 0.0,-0.55, -0.1])/3600




BLR = np.array([0.21791792, 0.33963964, 0.32822823, 0.3966967 ])*2.2
BLR[-1] = 0.68   
BLR = BLR/2

PSF=[[0.48, 0.39],
     [0.54, 0.4],
     [0.61, 0.45]]

#PSF = PSF/2

limx={'2QZJ': [2.15, 2.35, 2.2361, 616]} #
limx['LBQS'] = [2.17,2.24, 2.203, 250]
limx['HB89'] = [2.2,2.3, 2.2560, 500]

itere= np.array([1])
for i in itere:    
    #ID = New['XID'][i]
            
    
    
    ID = Sample['ID'][i].strip()
    z= Sample['z'][i]
    Hal_band='K_sin'
        
    H_file = Sample['Halpha cube'][i].strip() 
    
    if ID=='HB89':
        IDp='HB8903'
    if ID=='2QZJ':
        IDp='2QZJ00'       
    if ID=='LBQS':
        IDp= 'LBQS01'
        
        
        
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
    

    size = Hal_size[ID]
        
    rds = float(size[0])
    fl = size[1]
    print ('Radius and extraction region ', rds, fl)
        
    storage_H = IFU.choose_pixels(storage_H, plot_it, rad= rds , flg = fl)
        
    storage_H = IFU.astrometry_correction_GAIA(storage_H)
        
        
        
    storage_H = IFU.stack_sky(storage_H, Hal_band, 0, expand=0)
        
    storage_H = IFU.D1_spectra_collapse(storage_H,z, 'H', plot_it, addsave='_tot')
        
    storage_H = IFU.fitting_collapse_Halpha(storage_H,z, 1, broad=1)
    
    
    cont=1
    broad=1
     
    
    storage = storage_H
    
# =============================================================================
#   IMAGES Simulatnaous
# =============================================================================
    
    f = plt.figure(figsize=(15,15))
    
    
    
    Image_narrow = pyfits.getdata(ph.MyPATH+'Results_storage/Halpha/'+ID+'_Nearest_spaxel_fit.fits')
    IFU_header = pyfits.getheader(ph.MyPATH+'Results_storage/Halpha/'+ID+'_Nearest_spaxel_fit.fits')   
    
    Hal_map = smooth(Image_narrow[0],1.)/0.01
    
    use = np.isnan(Image_narrow[0])
    Hal_map[use] = np.nan
        
    shapes = np.shape(Hal_map)
    
        
    IFU_wcs= wcs.WCS(IFU_header).celestial
    
    ax1 = f.add_axes([0.05, 0.04, 0.22,0.22], projection= IFU_wcs)

    
    cmap = plt.get_cmap('cividis')
    cmap.set_bad(color='black')
    
    
    if ID=='2QZJ':
        
        Hal_map = Hal_map*1e16
    
    else:
        Hal_map = Hal_map*1e15
        print ('Normal normalisation')
    
    
    if ID=='2QZJ':
        fl = ax1.imshow(Hal_map, cmap=cmap, origin='low',vmax=7., aspect='auto', interpolation='nearest')
    
    
    else:
        fl = ax1.imshow(Hal_map, cmap=cmap, origin='low', aspect='auto', interpolation='nearest')
    
    
    spax = storage['Signal_mask'].copy()
    spax = spax[0,:,:]
    
    mask_nan =  np.ma.masked_invalid(Hal_map).mask   
    spax_tot = np.logical_or( mask_nan, spax)
    
       
    masked_flux = np.ma.array(data= Hal_map, mask= spax_tot)
    min_flux = np.ma.min(masked_flux)
    max_flux = np.ma.max(masked_flux)
    
    
    axcbar0 = plt.axes([ 0.27, 0.04 , 0.02, 0.22]) #plt.axes([0.055+ 0.2469,0.397,0.189,0.03])
    axcbar0.tick_params(direction='in')        
    cbar0 = f.colorbar(fl, cax=axcbar0 ,orientation='vertical')#, ticks= [min_flux,(min_flux+max_flux)/2, max_flux])
    #axcbar0.tick_params(axis='y',left='off',labelleft='off',right='on',labelright='on')
    #axcbar0.tick_params(axis='x',bottom='off',labelbottom='off',top='off',labeltop='off')
    
    if ID=='2QZJ':
        
        axcbar0.set_ylabel(r'SB (10$^{-16}$ erg/s/cm$^{2}$/arcsec$^{2}$)')
    
    else:
        axcbar0.set_ylabel(r'SB (10$^{-15}$ erg/s/cm$^{2}$/arcsec$^{2}$)')
        
    
    #axcbar0.xaxis.set_label_position('top')
    
    #axcbar0.tick_params(axis='y', direction='in', color='white')
    
    
    coor = storage['Median_stack_white_Center_data'][1:3].copy()
    #coor = np.array(coor, dtype=int)
    lim = lims[ID]
        
    ax1.set_xlim(coor[0]+lim[0],coor[0]+lim[1])
    ax1.set_ylim(coor[1]+lim[2],coor[1]+lim[3])
        
    ax1.set_autoscale_on(False)
    
    arc = np.round(1./(IFU_header['CDELT2']*3600))
    step = int(arc/4)   
    cent = np.round(storage['Median_stack_white_Center_data'][1:3]).copy()
    cent = np.array(cent, dtype=int)
    cent = cent-0.5
    
    alp=0.4
    ax1.plot( (cent[0]-step, cent[0]-step), (cent[1]+3*step, cent[1]-3*step), 'w-', alpha=alp)
    ax1.plot( (cent[0]+step, cent[0]+step), (cent[1]+3*step, cent[1]-3*step), 'w-', alpha=alp)           
    ax1.plot( (cent[0]-3*step, cent[0]-3*step), (cent[1]+3*step, cent[1]-3*step), 'w-', alpha=alp)
    ax1.plot( (cent[0]+3*step, cent[0]+3*step), (cent[1]+3*step, cent[1]-3*step), 'w-', alpha=alp)
       
    ax1.plot( (cent[0]-3*step, cent[0]+3*step), (cent[1]-3*step, cent[1]-3*step), 'w-', alpha=alp)
    ax1.plot( (cent[0]-3*step, cent[0]+3*step), (cent[1]+1*step, cent[1]+1*step), 'w-', alpha=alp)     
    ax1.plot( (cent[0]-3*step, cent[0]+3*step), (cent[1]-1*step, cent[1]-1*step), 'w-', alpha=alp)
    ax1.plot( (cent[0]-3*step, cent[0]+3*step), (cent[1]+3*step, cent[1]+3*step), 'w-', alpha=alp)
    
    
    
    #emplot.new_ALMA_contour_plot(storage, ax, both=0, prj = '4QSO'   )
    
    ax1.plot(coor[0], coor[1], 'r*', markersize=10)    
    
    dx = lim[1] - lim[0]
    dy = lim[3] - lim[2]
    x = coor[0]+lim[0]
    y = coor[1]+lim[2]
        
    ax1.text(dx*0.05+x, dy*0.9+y, IDp, color='white')
    
    deg_per_pix = IFU_header['CDELT2']
    arc_per_pix = deg_per_pix*3600
    
    lim_sc = lim*arc_per_pix
    
    from matplotlib.patches import Circle, Ellipse
    
    
    c = Circle((coor[0]+lim[1]*0.7, coor[1]-lim[3]*0.3), BLR[i]/arc_per_pix, edgecolor='limegreen', facecolor='limegreen', alpha=0.5)  
    ax1.add_patch(c)
    
    PSFs = PSF[0]
    
    c = Ellipse((coor[0]+lim[1]*0.7, coor[1]+lim[3]*0.7), PSFs[0]/arc_per_pix,PSFs[1]/arc_per_pix, 0, edgecolor='firebrick', facecolor='firebrick', alpha=0.5)  
    ax1.add_patch(c)

    emplot.new_ALMA_contour_plot(storage, ax1, both=0, prj = '4QSO'   )
    emplot.overide_axes_labels(f,ax1,(lim_sc), labelx=1, labely=1,tickin=1)
    
# =============================================================================
#  QSO sub    
# =============================================================================        
    
    Image_narrow = np.loadtxt(ph.MyPATH+'Results_storage/Sub_qso_map_hal/'+ID+'.txt')
    IFU_header = pyfits.getheader(ph.MyPATH+'Results_storage/Halpha/'+ID+'_Nearest_spaxel_fit.fits')   
        
    Hal_map = smooth(Image_narrow,1.)/0.01
    #Hal_map = Image_narrow
    IFU_wcs= wcs.WCS(IFU_header).celestial
    
   
    ax = f.add_axes([0.38, 0.04, 0.22,0.22], projection= IFU_wcs)
    
    #cmap = plt.get
    #cmap.set_bad(color='black')
    #cmap.set_under('k')
    
    
    
    
    if ID=='2QZJ':
        fl = ax.imshow(Hal_map*1e16, cmap=cmap, origin='low',vmin=0.5 , vmax=6., aspect='auto', interpolation='nearest')
        
    else:
        fl = ax.imshow(Hal_map*1e15, cmap=cmap, origin='low',vmin=0.2, aspect='auto', interpolation='nearest')
    
    spax = storage['Signal_mask'].copy()
    spax = spax[0,:,:]
    
    mask_nan =  np.ma.masked_invalid(Hal_map).mask   
    spax_tot = np.logical_or( mask_nan, spax)
    
       
    masked_flux = np.ma.array(data= Hal_map, mask= spax_tot)
    min_flux = np.ma.min(masked_flux)
    max_flux = np.ma.max(masked_flux)
    
    
    axcbar0 = plt.axes([ 0.6, 0.04 , 0.02, 0.22])  #plt.axes([0.055+ 0.2469,0.397,0.189,0.03])
    axcbar0.tick_params(direction='in')        
    cbar0 = f.colorbar(fl, cax=axcbar0 ,orientation='vertical')#, ticks= [min_flux,(min_flux+max_flux)/2, max_flux])
    #axcbar0.tick_params(axis='y',left='off',labelleft='off',right='on',labelright='on')
    #axcbar0.tick_params(axis='x',bottom='off',labelbottom='off',top='off',labeltop='off')
    
    if ID=='2QZJ':
        
        axcbar0.set_ylabel(r'SB (10$^{-16}$ erg/s/cm$^{2}$/arcsec$^{2}$)')
    
    else:
        axcbar0.set_ylabel(r'SB (10$^{-15}$ erg/s/cm$^{2}$/arcsec$^{2}$)')
        
    
    axcbar0.xaxis.set_label_position('top')
    
    axcbar0.tick_params(axis='x', direction='in', color='white')
    
    
    coor = storage['Median_stack_white_Center_data'][1:3].copy()
    #coor = np.array(coor, dtype=int)
    lim = lims[ID]
        
    ax.set_xlim(coor[0]+lim[0],coor[0]+lim[1])
    ax.set_ylim(coor[1]+lim[2],coor[1]+lim[3])
        
    ax.set_autoscale_on(False)
    
    arc = np.round(1./(IFU_header['CDELT2']*3600))
    step = int(arc/4)   
    cent = np.round(storage['Median_stack_white_Center_data'][1:3]).copy()
    cent = np.array(cent, dtype=int)
    cent = cent-0.5
    
       
    alp=0.4
    ax.plot( (cent[0]-step, cent[0]-step), (cent[1]+3*step, cent[1]-3*step), 'w-', alpha=alp)
    ax.plot( (cent[0]+step, cent[0]+step), (cent[1]+3*step, cent[1]-3*step), 'w-', alpha=alp)           
    ax.plot( (cent[0]-3*step, cent[0]-3*step), (cent[1]+3*step, cent[1]-3*step), 'w-', alpha=alp)
    ax.plot( (cent[0]+3*step, cent[0]+3*step), (cent[1]+3*step, cent[1]-3*step), 'w-', alpha=alp)
       
    ax.plot( (cent[0]-3*step, cent[0]+3*step), (cent[1]-3*step, cent[1]-3*step), 'w-', alpha=alp)
    ax.plot( (cent[0]-3*step, cent[0]+3*step), (cent[1]+1*step, cent[1]+1*step), 'w-', alpha=alp)     
    ax.plot( (cent[0]-3*step, cent[0]+3*step), (cent[1]-1*step, cent[1]-1*step), 'w-', alpha=alp)
    ax.plot( (cent[0]-3*step, cent[0]+3*step), (cent[1]+3*step, cent[1]+3*step), 'w-', alpha=alp)
    
    
    
    #emplot.new_ALMA_contour_plot(storage, ax, both=0, prj = '4QSO'   )
    
    ax.plot(coor[0], coor[1]-0.5, 'r*', markersize=10)    
    
    dx = lim[1] - lim[0]
    dy = lim[3] - lim[2]
    x = coor[0]+lim[0]
    y = coor[1]+lim[2]
        
    ax.text(dx*0.05+x, dy*0.9+y, IDp, color='white')
    
    deg_per_pix = IFU_header['CDELT2']
    arc_per_pix = deg_per_pix*3600
    
    lim_sc = lim*arc_per_pix
    
    c = Circle((coor[0]+lim[1]*0.7, coor[1]-lim[3]*0.3), BLR[i]/arc_per_pix, edgecolor='limegreen', facecolor='limegreen', alpha=0.5)  
    ax1.add_patch(c)
    
    PSFs = PSF[0]
    
    c = Ellipse((coor[0]+lim[1]*0.7, coor[1]+lim[3]*0.7), PSFs[0]/arc_per_pix,PSFs[1]/arc_per_pix, 0, edgecolor='firebrick', facecolor='firebrick', alpha=0.5)  
    ax1.add_patch(c)
    
    emplot.new_ALMA_contour_plot(storage, ax, both=0, prj = '4QSO'   )
    
    emplot.overide_axes_labels(f,ax,(lim_sc), labelx=1, labely=1,tickin=1)
    
    
# =============================================================================
#     Residual 
# =============================================================================
    Data = pyfits.getdata(PATH+'Four_Quasars/Results_storage/Hal_fitting/'+ID+'_Nearest_spaxel_fit_res.fits')
    header = pyfits.getheader(PATH+'Four_Quasars/Results_storage/Hal_fitting/'+ID+'_Nearest_spaxel_fit_res.fits')
    
    axres = f.add_axes([0.71, 0.04, 0.22,0.22], projection= IFU_wcs)
    
    outo = storage_H['1D_fit_Halpha_mul']
    wave= storage['obs_wave']
    
    halc = outo.params['Han_center'].value
    rng = 500./halc*3e5#limx[ID][3]/3e5*limx[ID][2]
    
    #collapse = np.where((wave<limx[ID][2]+rng/2) & (wave>limx[ID][2]-rng/2))[0]
    
    collapse = np.where((wave<halc+rng/2) & (wave>halc-rng/2))[0]
     
    Res_image = Data[collapse,:,:].sum(axis=0)
    
    #Res_image = smooth(Res_image, 1)
    
    Res_image[np.where(Res_image==0)] = np.nan
      
    rms = np.nanstd(Res_image- np.nanmean(Res_image))
    
    
    fl = axres.imshow(Res_image/rms,vmin=-3, vmax=3,cmap=cmap,   origin='low')
    
    axres.contour(Res_image/rms, levels=(3, 5), linestyles='solid', colors='red')
    
    
    im = lims[ID]
    
    axres.set_xlim(coor[0]+lim[0],coor[0]+lim[1])
    axres.set_ylim(coor[1]+lim[2],coor[1]+lim[3])
    
    deg_per_pix = header['CDELT2']
    arc_per_pix = deg_per_pix*3600
    
    dx = lim[1] - lim[0]
    dy = lim[3] - lim[2]
    x = coor[0]+lim[0]
    y = coor[1]+lim[2]
        
    
    lim_sc = lim*arc_per_pix
    
    emplot.overide_axes_labels(f,axres,(lim_sc), labelx=1, labely=1,tickin=1)
    
    axres.text(dx*0.05+x, dy*0.9+y, IDp, color='white')
    
    axcbar0 = plt.axes([ 0.92, 0.04 , 0.02, 0.22])  #plt.axes([0.055+ 0.2469,0.397,0.189,0.03])
    axcbar0.tick_params(direction='in')        
    cbar0 = f.colorbar(fl, cax=axcbar0 ,orientation='vertical')#, 
    
    axcbar0.set_ylabel(r'Residual ($\sigma$)')
    
    
    axres.plot( (cent[0]-step, cent[0]-step), (cent[1]+3*step, cent[1]-3*step), 'w-', alpha=alp)
    axres.plot( (cent[0]+step, cent[0]+step), (cent[1]+3*step, cent[1]-3*step), 'w-', alpha=alp)           
    axres.plot( (cent[0]-3*step, cent[0]-3*step), (cent[1]+3*step, cent[1]-3*step), 'w-', alpha=alp)
    axres.plot( (cent[0]+3*step, cent[0]+3*step), (cent[1]+3*step, cent[1]-3*step), 'w-', alpha=alp)
       
    axres.plot( (cent[0]-3*step, cent[0]+3*step), (cent[1]-3*step, cent[1]-3*step), 'w-', alpha=alp)
    axres.plot( (cent[0]-3*step, cent[0]+3*step), (cent[1]+1*step, cent[1]+1*step), 'w-', alpha=alp)     
    axres.plot( (cent[0]-3*step, cent[0]+3*step), (cent[1]-1*step, cent[1]-1*step), 'w-', alpha=alp)
    axres.plot( (cent[0]-3*step, cent[0]+3*step), (cent[1]+3*step, cent[1]+3*step), 'w-', alpha=alp)
    
# =============================================================================
#   SPECTRUMS
# =============================================================================
    
    center =  storage['Median_stack_white_Center_data'][1:3].copy()
    shapes = storage['dim']
    ID = storage['X-ray ID']
    wave= storage['obs_wave']
    z = storage['z_guess']    
    
    
    # Creating a mask for all spaxels. 
    mask_catch = storage['flux'].mask.copy()
    
    mask = mask_catch.copy()
    
    header  = storage['header']
    #arc = np.round(1./(header['CD2_2']*3600))
    arc = np.round(1./(header['CDELT2']*3600))
    
    step = int(arc/4)
    
    cent = np.round(storage['Median_stack_white_Center_data'][1:3])
    cent = np.array(cent, dtype=int)
    
    
   
    indxs = np.array([1,0,-1])
    indys = np.array([-1,0,1])
    
    
    indx=0
    indy=0
    
    
    mask[:,:,:] = True          
    mask[:,cent[1]+indx*step*2- step:cent[1]+indx*step*2+ step,cent[0]+indy*step*2- step:cent[0]+indy*step*2+ step]= False
          
            
    flux = np.ma.array(data=storage['flux'].data, mask= mask) 
    
    Spectra = np.ma.sum(flux, axis=(1,2))
    Spectra = np.ma.array(data=Spectra.data, mask=storage['sky_clipped_1D'])
    error = IFU.STD_calc(wave/(1+z)*1e4, Spectra, 'H')* np.ones_like(Spectra)
     
    
    Save_spec = np.zeros((4,len(Spectra)))
    
    Save_spec[0,:] = wave
    Save_spec[1,:] = Spectra
    Save_spec[2,:] = error
    Save_spec[3,:] = storage['sky_clipped_1D']
    
    
    if ID=='HB89':
        outo = storage_H['1D_fit_Halpha_mul']
    else:
        outo, chi2 = emfit.fitting_Halpha_mul_testing(wave,Spectra,error,z, broad=1)
        
    #
    cmap = plt.cm.viridis
    
    Ha_flux = np.array([IFU.flux_measure_ind(outo, wave, 'H', use='BPT')])
    
    np.savetxt(ph.MyPATH+'Results_storage/'+ID+'_Halpha_nuclear_flux.txt',Ha_flux)
    
    for ix in np.array([0,1,2]):
        for iy in np.array([0,1,2]):
            
            
            indx = indxs[ix]
            indy = indys[iy]
            
            print( indx, indy)
            
            mask[:,:,:] = True
            
            mask[:,cent[1]+indx*step*2- step:cent[1]+indx*step*2+ step,cent[0]+indy*step*2- step:cent[0]+indy*step*2+ step]= False
            
            
            ax = f.add_axes([0.37+ indy*0.32, 0.56+ indx*0.23, 0.29,0.2])
            
            
            #axes_im[ix,iy].imshow(storage['Median_stack_white'], origin='low', vmax=0.211, alpha=0.2)    
            
            #axes_im[ix,iy].imshow(np.ma.array(data=storage['Median_stack_white'], mask=mask[0,:,:]), origin='low', vmax=0.211)    
           
            
            flux = np.ma.array(data=storage['flux'].data, mask= mask) 
    
            Spectra = np.ma.sum(flux, axis=(1,2))
            Spectra = np.ma.array(data=Spectra.data, mask=storage['sky_clipped_1D'])
            error = IFU.STD_calc(wave/(1+z)*1e4, Spectra, 'H')* np.ones_like(Spectra)
            
            if (ix==1) & (iy==1) & (ID=='HB89'):
                
                Spectra = storage['1D_spectrum'].copy()
                error = storage['1D_spectrum_er'].copy()
                print('Overiride')
    
            
            if ID=='2QZJ':
                if (ix==1) & (iy==2):
                    tt=0.08
                else:
                    tt='n'
                
                if (ix==1) & (iy==1):
                    out, chi2 = emfit.fitting_Halpha_mul_testing(wave,Spectra,error,z, broad=broad, cont=cont, decompose=outo) 
                else:
                    out, chi2 = emfit.fitting_Halpha_mul_bkp_2QZJ(wave,Spectra,error,z, broad=broad, cont=cont, decompose=outo, cont_norm=tt) 
                
            
            
            else:
                out, chi2 = emfit.fitting_Halpha_mul_testing(wave,Spectra, error,z, broad=broad, cont=cont, decompose=outo)
                
            if (ix==1) & (iy==1) & (ID=='HB89'):
                out=outo
            
    
            if ID=='2QZJ':
                
                if (indx==0) :
                    emplot.plotting_Halpha(wave, Spectra, ax, out, 'mul',z, title=0, yscale='model')
                    
                    
                    SNR, chi2 = IFU.SNR_calc(Spectra, wave, out, 'H',z, mul=1)
            
                    hal_fl  = IFU.flux_measure_ind(out,wave, 'H', use='BPT', error=0)
            
                    print (hal_fl)
                    print ('Hal/NII ', np.log10(hal_fl/IFU.flux_measure_ind(out,wave, 'NII', error=0)))
                    
                    
                    if indy==0:
                        ymin = -0.1
                        ymax = 0.7
                    
                    else:
                        ymin, ymax = ax.get_ylim()
                        ax.set_ylim(1.5*ymin, ymax*1.5)
                    
                    
                    ymin, ymax = ax.get_ylim()
                    ax.fill_between((6710, 6736), y1 =-0.1*max(out.eval(x=wave)), y2 = 100 , color='grey', alpha=0.4)  
                    #ax.text(6737, ymax*0.5, '[SII]', color='firebrick')
                       
                    ax.text(6682, ymax*0.7  , 'Narrow H$\\alpha$ flux: \n%.1f $x 10^{-16}$ erg/s/cm$^{2}$' %(hal_fl*1e16))
                    ax.set_xlim(6200,7000)
                    
                    # Saving the data to use elsewhere 
                    Save_spec = {'wave' : wave}
                    Save_spec['data'] = Spectra
                    Save_spec['error'] = error
                    Save_spec['z'] = z
                    Save_spec['total'] = out.eval(x=wave)
                    
                    Save_spec['cont'] = out.eval_components(x=wave)['linear']
                    
                    Save_spec['Haw'] = out.eval_components(x=wave)['Haw_']        
                    Save_spec['Han'] = out.eval_components(x=wave)['Han_']
                    
                    Save_spec['Nr'] = out.eval_components(x=wave)['Nr_']
                    Save_spec['Nb'] = out.eval_components(x=wave)['Nb_']
                    try:
                        Save_spec['Hao'] = out.eval_components(x=wave)['Hao_']
                    except:
                        ldsf=999999
                    
                    Save_spec['Hal_flx'] = hal_fl
                    
                    vels = (wave-outo.params['Han_center'].value)/outo.params['Han_center'].value*3e5
                
                    axres = f.add_axes([0.39+ indy*0.32+0.005, 0.56+0.13+ indx*0.23, 0.29/4,0.2/3])
                    axres.plot(vels, Spectra.data-out.eval(x=wave), color='grey', drawstyle='steps-mid', alpha=0.2)
                    axres.set_ylabel('Flux (rms )', fontsize=9)
                    
                    flux_sc = Spectra.data[np.invert(Spectra.mask)]
                    wave_sc= wave[np.invert(Spectra.mask)]
            
                    vels_sc = (wave_sc-outo.params['Han_center'].value)/outo.params['Han_center'].value*3e5
                    
                    use = np.where((vels_sc<2000) & (vels_sc>-2000))
                    
                    er = np.std(flux_sc-out.eval(x=wave_sc))
                    axres.plot(vels_sc,flux_sc-out.eval(x=wave_sc), drawstyle='steps-mid')
                    
                    axres.set_xlim(-990,990)
                    try:     
                        axres.set_ylim(-4*er, 4*er)
                    except:
                        er = 0.18
                        axres.set_ylim(-4*er, 4*er)
                        
                    axres.set_yticks([0,3*er])
                    axres.set_xticks([-700,0,700])
                    
                    labels = [item.get_text() for item in axres.get_yticklabels()]
                    labels[1] =  3
                    labels[0] =  0
                    
                    axres.set_yticklabels(labels)
                    axres.tick_params(direction='in')
                        
                    
                    emplot.save_obj(Save_spec, 'Results_storage/Spectrums/Regions/'+ID+'_Hal_'+str(iy)+'_'+str(ix)+'_mod')
                
                elif (indx==1) & (indy==0):
                    emplot.plotting_Halpha(wave, Spectra, ax, out, 'mul',z, title=0, yscale='model')
                    
                    
                    SNR, chi2 = IFU.SNR_calc(Spectra, wave, out, 'H',z, mul=1)
            
                    hal_fl  = IFU.flux_measure_ind(out,wave, 'H', use='BPT', error=0)
            
                    print (hal_fl)
                    
                    print ('Hal/NII ',np.log10(hal_fl/IFU.flux_measure_ind(out,wave, 'NII', error=0)))
                    ymin, ymax = ax.get_ylim()
                    ax.set_ylim(1.5*ymin, ymax*1.5)
                    ymin, ymax = ax.get_ylim()
                    
                    ax.fill_between((6710, 6736), y1 =-0.1*max(out.eval(x=wave)), y2 =100 , color='grey', alpha=0.4)  
                    #ax.text(6737, ymax*0.5, '[SII]', color='firebrick')
                       
                    ax.text(6682, ymax*0.7  , 'Narrow H$\\alpha$ flux: \n%.1f $x 10^{-16}$ erg/s/cm$^{2}$' %(hal_fl*1e16))
                    ax.set_xlim(6200,7000)
                    
                    # Saving the data to use elsewhere 
                    Save_spec = {'wave' : wave}
                    Save_spec['data'] = Spectra
                    Save_spec['error'] = error
                    Save_spec['z'] = z
                    Save_spec['total'] = out.eval(x=wave)
                    
                    Save_spec['cont'] = out.eval_components(x=wave)['linear']
                    
                    Save_spec['Haw'] = out.eval_components(x=wave)['Haw_']        
                    Save_spec['Han'] = out.eval_components(x=wave)['Han_']
                    
                    Save_spec['Nr'] = out.eval_components(x=wave)['Nr_']
                    Save_spec['Nb'] = out.eval_components(x=wave)['Nb_']
                    
                    Save_spec['Hal_flx'] = hal_fl
                    
                    
                    vels = (wave-outo.params['Han_center'].value)/outo.params['Han_center'].value*3e5
                
                
                    axres.plot(vels, Spectra.data-out.eval(x=wave), color='grey', drawstyle='steps-mid', alpha=0.2)
                    
                    
                    flux_sc = Spectra.data[np.invert(Spectra.mask)]
                    wave_sc= wave[np.invert(Spectra.mask)]
            
                    vels_sc = (wave_sc-outo.params['Han_center'].value)/outo.params['Han_center'].value*3e5
                    
                    use = np.where((vels_sc<2000) & (vels_sc>-2000))
                    
                    vels = (wave-outo.params['Han_center'].value)/outo.params['Han_center'].value*3e5
                
                    axres = f.add_axes([0.39+ indy*0.32+0.005, 0.56+0.13+ indx*0.23, 0.29/4,0.2/3])
                    axres.plot(vels, Spectra.data-out.eval(x=wave), color='grey', drawstyle='steps-mid', alpha=0.2)
                    axres.set_ylabel('Flux (rms )', fontsize=9)
                    
                    flux_sc = Spectra.data[np.invert(Spectra.mask)]
                    wave_sc= wave[np.invert(Spectra.mask)]
            
                    vels_sc = (wave_sc-outo.params['Han_center'].value)/outo.params['Han_center'].value*3e5
                    
                    use = np.where((vels_sc<2000) & (vels_sc>-2000))
                    
                    er = np.std(flux_sc-out.eval(x=wave_sc))
                    
        
                    axres.plot(vels_sc,flux_sc-out.eval(x=wave_sc), drawstyle='steps-mid')
                    
                    try:     
                        axres.set_ylim(-2*er, 2*er)
                    except:
                        er = 0.18
                        axres.set_ylim(-2*er, 2*er)
                        
                    axres.set_yticks([0,1.5*er])
                    axres.set_xticks([-700,0,700])
                    axres.set_xlim(-990,990)
                    
                    labels = [item.get_text() for item in axres.get_yticklabels()]
                    labels[1] =  3
                    labels[0] =  0
                    
                    axres.set_yticklabels(labels)
                    axres.tick_params(direction='in')
                    
                    emplot.save_obj(Save_spec, 'Results_storage/Spectrums/Regions/'+ID+'_Hal_'+str(iy)+'_'+str(ix)+'_mod')
                
                else:
                    wv_rest = wave/(1+z)*1e4
                    fit_loc = np.where((wv_rest>6000.)&(wv_rest<7500.))[0]
                    
                    ax.plot(wv_rest[fit_loc], Spectra.data[fit_loc], color='grey', drawstyle='steps-mid', alpha=0.2)
   
                    
                    
                    flux = Spectra.data[np.invert(Spectra.mask)]
                    wv_rst_sc= wv_rest[np.invert(Spectra.mask)]
        
                    fit_loc_sc = np.where((wv_rst_sc>6000.)&(wv_rst_sc<6999.))[0]   
                    
                    #from scipy.signal import medfilt as mdf
                    ax.plot(wv_rst_sc[fit_loc_sc],flux[fit_loc_sc], drawstyle='steps-mid')
                    
                    ax.set_ylim(-0.5,0.5)
                
                    vels = (wave-outo.params['Han_center'].value)/outo.params['Han_center'].value*3e5
                
                
                    axres.plot(vels, Spectra.data-out.eval(x=wave), color='grey', drawstyle='steps-mid', alpha=0.2)
                    
                    
                    flux_sc = Spectra.data[np.invert(Spectra.mask)]
                    wave_sc= wave[np.invert(Spectra.mask)]
            
                    vels_sc = (wave_sc-outo.params['Han_center'].value)/outo.params['Han_center'].value*3e5
                    
                    use = np.where((vels_sc<2000) & (vels_sc>-2000))
                    
                    
        
                    
            else:
                emplot.plotting_Halpha(wave, Spectra, ax, out, 'mul',z, title=0)
                
                   
                SNR, chi2 = IFU.SNR_calc(Spectra, wave, out, 'H',z, mul=1)
            
                hal_fl  = IFU.flux_measure_ind(out,wave, 'H', use='BPT', error=0)
                
                ax.fill_between((6710, 6736), y1 =-0.1*max(out.eval(x=wave)), y2 = 1.2*max(out.eval(x=wave)) , color='grey', alpha=0.4)  
                #ax.text(6737, ymax*0.5, '[SII]', color='firebrick')
            
                print (hal_fl)
                print ('Hal/NII ',np.log10(hal_fl/IFU.flux_measure_ind(out,wave, 'NII', error=0)))
                ymin, ymax = ax.get_ylim()
                       
                ax.text(6620, ymax*0.7  , 'Narrow H$\\alpha$ flux: \n%.1f $x 10^{-16}$ erg/s/cm$^{2}$' %(hal_fl*1e16))
                ax.set_xlim(6200,7000)
                
                # Saving the data to use elsewhere 
                Save_spec = {'wave' : wave}
                Save_spec['data'] = Spectra
                Save_spec['error'] = error
                Save_spec['z'] = z
                Save_spec['total'] = out.eval(x=wave)
                
                Save_spec['cont'] = out.eval_components(x=wave)['linear']
                
                Save_spec['Haw'] = out.eval_components(x=wave)['Haw_']        
                Save_spec['Han'] = out.eval_components(x=wave)['Han_']
                
                Save_spec['Nr'] = out.eval_components(x=wave)['Nr_']
                Save_spec['Nb'] = out.eval_components(x=wave)['Nb_']
                
                try:
                    Save_spec['Hao'] = out.eval_components(x=wave)['Hao_']    
                except:
                    lds=9999
                
                
                vels = (wave-outo.params['Han_center'].value)/outo.params['Han_center'].value*3e5
                
                axres = f.add_axes([0.39+ indy*0.32+0.005, 0.56+0.13+ indx*0.23, 0.29/4,0.2/3])
                axres.plot(vels, Spectra.data-out.eval(x=wave), color='grey', drawstyle='steps-mid', alpha=0.2)
                axres.set_ylabel('Flux (rms )', fontsize=9)
                
                
                flux_sc = Spectra.data[np.invert(Spectra.mask)]
                wave_sc= wave[np.invert(Spectra.mask)]
        
                vels_sc = (wave_sc-outo.params['Han_center'].value)/outo.params['Han_center'].value*3e5
                
                use = np.where((vels_sc<2000) & (vels_sc>-2000))
                
                er = np.std(flux_sc-out.eval(x=wave_sc))
                
    
                axres.plot(vels_sc,flux_sc-out.eval(x=wave_sc), drawstyle='steps-mid')
                
                axres.hlines(3*er, -3000,3000, color='k', linestyle='dashed')
                
                axres.text(-750, 2*er,r'3$\sigma$ limit' )
                
                
                axres.set_xlim(-990,990)
                axres.set_ylim(-4*er, 4*er)
                
                axres.set_yticks([0,3*er])
                axres.set_xticks([-700,0,700])
                
                labels = [item.get_text() for item in axres.get_yticklabels()]
                labels[1] =  3
                labels[0] =  0
                
                axres.set_yticklabels(labels)
                axres.tick_params(direction='in')
                
                Save_spec['Hal_flx'] = hal_fl
                
                emplot.save_obj(Save_spec, 'Results_storage/Spectrums/Regions/'+ID+'_Hal_'+str(iy)+'_'+str(ix)+'_mod')
                   
                    
            ax.set_xlim(6300,6900)
            #ax.text(6682, ymax*0.1  , str(int(cent[0]+indy*step*2))+', '+ str(int(cent[1]+indx*step*2)) )
            if indx==-1:
                ax.set_xlabel('Rest-frame wavelength (ang)', fontsize=12)
                #axres.set_xlabel('Velocity (km/s)', fontsize=12)
            
            
            if indy==-1:
                ax.set_ylabel(r'x$10^{-17}$ (erg/s/cm$^2$/ang)', fontsize=12)
                #axres.set_ylabel(r'x$10^{-17}$ (erg/s/cm$^2$/ang)', fontsize=12)
    
    cent = cent-0.5
    
    

    
    
    #f.savefig(ph.MyPATH+'Graphs/Paper_plots/Fig9-11_'+ID+'_Halpha_regions_resinset.pdf')#, bbox_inches = 'tight')
    

plt.show()
        
        
     
