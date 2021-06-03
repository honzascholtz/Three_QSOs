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


OIII= False
Hal = True
Hbeta = False

Spat = False


def smooth(image,sm):
    
    from astropy.convolution import Gaussian2DKernel
    from scipy.signal import convolve as scipy_convolve
    from astropy.convolution import convolve
    
    gauss_kernel = Gaussian2DKernel(sm)

    con_im = convolve(image, gauss_kernel)
    
    #con_im = con_im#*image/image
    
    return con_im  
 
###############################################################################
# Halpha nad ALMA
###############################################################################
column = 0
row = 0

center = {'LBQS' : np.array([53,47])}
center['HB89'] = np.array([53.7,47.9])
center['2QZJ'] = np.array([53.0, 47.5])

lims={'HB89':  np.array([-10,10,-10,10])}
lims['2QZJ'] = np.array([-10,12,-10,12])
lims['LBQS'] = np.array([-10,10,-10,10])


OIII_size={'HB89': np.array(['0.2', '1'], dtype=str)}
OIII_size['LBQS'] = np.array(['0.2', '1'], dtype=str)
OIII_size['2QZJ'] = np.array(['0.2', '1'], dtype=str)



IFU_psfs = {'HB89': np.array([-1.1,-0.6,-0.08])/3600}
IFU_psfs['2QZJ'] = np.array([-0.8,-0.7, -0.1 ])/3600
IFU_psfs['LBQS'] = np.array([ 0.0,-0.55, -0.1])/3600


itere=np.array([1])
for i in itere:    
    #ID = New['XID'][i]
              
    ID = Sample['ID'][i].strip()
    
    OIII_band='Hsin'
    
    H_file=ID+'_H_fl'

    ID = Sample['ID'][i].strip()
    z= Sample['z'][i]
    if ID=='2QZJ':
        z= 2.4063
        
    print (ID)
    
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

    storage_O= IFU.fitting_collapse_OIII(storage_O,z, 0)
    
    storage= storage_O
    
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
    print (step)
    
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
    error = IFU.STD_calc(wave/(1+z)*1e4, Spectra, 'OIII')* np.ones_like(Spectra.data)  
    
    hbw=4861
    if (ID=='HB89'):
        hbw=4870.
    if (ID=='LBQS'):
        hbw=4872
    outo = emfit.fitting_OIII_Hbeta_qso_mul(wave,Spectra, error,z, hbw=hbw)
    Saves = emplot.load_obj('Four_Quasars/Results_storage/Spectrums/Regions/'+ID+'_OIII_'+str(1)+'_'+str(1))
    
    wave = Saves['wave']#/(1+z)*1e4              
    fit = Saves['total']  
    
    wp = wave[np.argmax(fit)]
    z = (wp*1e4/5008)-1
    print(z)
    
    
    f = plt.figure(figsize=(15,10))
    
    #gd, axes_im = plt.subplots(3,3, figsize=(15,15))
    
    for ix in np.array([0,1,2]):
        for iy in np.array([0,1,2]):
            
            
            indx = indxs[ix]
            indy = indys[iy]
            print(' ')
            print(' ')
            print (indx, indy)
            print('Central coor', int(cent[0]+indy*step*2), int(cent[1]+indx*step*2))
            
            mask[:,:,:] = True
            
            mask[:,cent[1]+indx*step*2- step:cent[1]+indx*step*2+ step,cent[0]+indy*step*2- step:cent[0]+indy*step*2+ step]= False
            
            
            ax = f.add_axes([0.37+ indy*0.32, 0.56+ indx*0.23, 0.29,0.2])
            
# =============================================================================
#             Testing the region extraction
# =============================================================================
            
            
            ax.axvline(wp/(1+z)*1e4, linestyle='--', color='k', alpha=0.3)
            
            Saves = emplot.load_obj('Four_Quasars/Results_storage/Spectrums/Regions/'+ID+'_OIII_'+str(iy)+'_'+str(ix))
    
            wave = Saves['wave']/(1+z)*1e4
            Spec = Saves['data']
            
            fit = Saves['total']
            Hbw = Saves['Hbw']
            Hbn = Saves['Hbn']
            
            o3rn = Saves['o3rn']
            o3rw = Saves['o3rw']
            
            o3bn = Saves['o3bn']
            o3bw = Saves['o3bw']
            
            
            flux = Spec.data[np.invert(Spec.mask)]
            wv_rst_sc= wave[np.invert(Spec.mask)]
            
            
            ax.plot(wave, Spec.data, drawstyle='steps-mid', linestyle='solid', color='grey', alpha=0.2)
            
            ax.plot(wv_rst_sc, flux, drawstyle='steps-mid', linestyle='solid')
            
            ax.plot(wave, (o3rn+o3bn), 'green', linestyle='dashed')
            ax.plot(wave, (o3rw+o3bw), 'blue', linestyle='dashed')
            ax.plot(wave, fit, 'red', linestyle='dashed')
            
            ax.plot(wave, Hbw, 'orange', linestyle='dotted')
            ax.plot(wave, Hbn, 'orange', linestyle='dashed')
                       
            ymax = max(fit)
            
            cr =0
            if ID=='HB89':
                off = (wp-Saves['wave'][np.argmax(o3rw)])/wp*3e5
            
            else:
                off = Saves['offset']
                
            w80 = Saves['W80']
            F_broad = Saves['Broad flux']
            
            
            
            if indx==-1:
                ax.set_xlabel('Rest-frame wavelength (ang)', fontsize=12)
            
            
            if indy==-1:
                ax.set_ylabel(r'x$10^{-17}$ (erg/s/cm$^2$/ang)', fontsize=12)
            
            ax.set_ylim(-0.1*ymax, 1.2*ymax)
            
            ymin, ymax = ax.get_ylim()
            
            if ID=='2QZJ':
                if ix==2:
                    ymax = ymax*2
                    
                elif (ix==0)&(iy==0):
                    ymax = ymax*2
                
                else:
                    ymax = ymax*1.5
                    
                ax.set_ylim(ymin, ymax)
            
            ax.text(4720, ymax*0.85  , r'W80 %.0f km s$^{-1}$' %(w80[2]))
            ax.text(4720, ymax*0.75  , r'Offset %.0f km s$^{-1}$' %(-off))
            ax.text(4720, ymax*0.65  , r'Broad Flux %.1f $\times 10^{-16}$ ergs/s' %(F_broad*1e16))
            
            ax.set_xlim(4700,5100)
            
            
            
            if indx==-1:
                ax.set_xlabel('Rest-frame wavelength (ang)', fontsize=12)
            
            
            if indy==-1:
                ax.set_ylabel(r'x$10^{-17}$ (erg/s/cm$^2$/ang)', fontsize=12)
            
            
 
# =============================================================================
#     Setting the velocity boundries
# =============================================================================
    alp=0.6
    
    st='w-'
    lnw = 2
    #z= Sample['z'][i]
    if ID=='2QZJ':
        z= 2.4063


# =============================================================================
#  Plotting their definition of outflows
# =============================================================================


    Maps = pyfits.getdata(PATH+'Four_Quasars/Results_storage/'+ID+'_OIII_mapsind.fits')
    Results = pyfits.getdata(PATH+'Four_Quasars/Results_storage/OIII/'+ID+'_Individual_spaxel_fit.fits')  
    IFU_header = pyfits.getheader(PATH+'Four_Quasars/Results_storage/OIII/'+ID+'_Individual_spaxel_fit.fits')   
    IFU_wcs= wcs.WCS(IFU_header).celestial
    
    cts=coor = center
    ax = f.add_axes([0.01, 0.04, 0.36, 0.23], projection=IFU_wcs)
    
    cnt =  5008*(1+z)/1e4
    med = (Results[3]-cnt)/cnt*3e5
    
    if ID=='LBQS':
        vn = -200
        vx = 300
    
        lvls = (-150,-100,-50)
    
    if ID=='HB89':
        vn = -300  #-700
        vx = -50    #-300
    
        lvls = (-200, -125) #(-550, -500)
        
        if ID=='HB89':
            print ('HB89 Cent ', cnt)
        
        med = med
        med[55,56:63] = -300+300 #-300
        med[54,58:63] = -350+300 #-350
        med[53,60:63] = -400+300 #-400
        
        med[53,60:63] = -400+300 #-400
        med[52,60:63] = -400+300 #-400
        med[51,58:63] = -400+300 #-400
        med[50,60:63] = -400+300 #-400
        #med[50,58:63] = -400
    
    if ID=='2QZJ':    
        
        vn=-550
        vx= -400
        lvls = (-500,-450)
    
    smt = smooth(med,1.5)

    use = np.isnan(med)
    smt[use] = np.nan
   
    cmap = plt.cm.coolwarm
    cmap.set_bad(color='black')
    
    fl = ax.imshow(smt,vmin=vn, vmax=vx, origin='low', cmap=cmap)
    ax.contour(smt, levels=lvls, colors='white', linewidth=2, linestyle='dashed')
    
    axcbar0 = f.add_axes([ 0.264, 0.04 , 0.02, 0.23]) #plt.axes([0.055+ 0.2469,0.397,0.189,0.03])
    axcbar0.tick_params(direction='in')        
    cbar0 = f.colorbar(fl, cax=axcbar0 ,orientation='vertical')
    
    axcbar0.set_ylabel('Broad offset (km/s)')
    
    lim = lims[ID]
    
    deg_per_pix = IFU_header['CDELT2']
    arc_per_pix = deg_per_pix*3600
    
    lim_sc = lim*arc_per_pix
    
    ax.set_xlim(coor[0]+lim[0],coor[0]+lim[1])
    ax.set_ylim(coor[1]+lim[2],coor[1]+lim[3])
    
    ax.plot( cts[0], cts[1], 'r*', markersize=12)
    
    emplot.overide_axes_labels(f,ax,(lim_sc), labelx=1, labely=1, tickin=2)
    
    ax.plot( (cent[0]-step, cent[0]-step), (cent[1]+3*step, cent[1]-3*step), st, alpha=alp, linewidth=lnw)
    ax.plot( (cent[0]+step, cent[0]+step), (cent[1]+3*step, cent[1]-3*step), st, alpha=alp, linewidth=lnw)           
    ax.plot( (cent[0]-3*step, cent[0]-3*step), (cent[1]+3*step, cent[1]-3*step), st, alpha=alp, linewidth=lnw)
    ax.plot( (cent[0]+3*step, cent[0]+3*step), (cent[1]+3*step, cent[1]-3*step), st, alpha=alp, linewidth=lnw)
       
    ax.plot( (cent[0]-3*step, cent[0]+3*step), (cent[1]-3*step, cent[1]-3*step), st, alpha=alp, linewidth=lnw)
    ax.plot( (cent[0]-3*step, cent[0]+3*step), (cent[1]+1*step, cent[1]+1*step), st, alpha=alp, linewidth=lnw)     
    ax.plot( (cent[0]-3*step, cent[0]+3*step), (cent[1]-1*step, cent[1]-1*step), st, alpha=alp, linewidth=lnw)
    ax.plot( (cent[0]-3*step, cent[0]+3*step), (cent[1]+3*step, cent[1]+3*step), st, alpha=alp, linewidth=lnw)
    


# =============================================================================
#      Loading Maps and results   
# =============================================================================   
    Maps = pyfits.getdata(PATH+'Four_Quasars/Results_storage/'+ID+'_OIII_mapsind.fits')
    Results = pyfits.getdata(PATH+'Four_Quasars/Results_storage/OIII/'+ID+'_Individual_spaxel_fit.fits')  
    IFU_header = pyfits.getheader(PATH+'Four_Quasars/Results_storage/OIII/'+ID+'_Individual_spaxel_fit.fits')   
    IFU_wcs= wcs.WCS(IFU_header).celestial
        
    cts=coor = center
    deg_per_pix = IFU_header['CDELT2']
    arc_per_pix = deg_per_pix*3600
    
    lim = lims[ID]    
    lim_sc = lim*arc_per_pix
    
    IFU_wcs= wcs.WCS(IFU_header).celestial


    
# =============================================================================
#     W80 Plot
# =============================================================================
    ax = f.add_axes([0.33, 0.04, 0.35, 0.23], projection=IFU_wcs)
    
    smt = smooth(Maps[2,:,:],1.5)
    
    use = np.isnan(Maps[2,:,:])
    smt[use] = np.nan
    
    if ID=='HB89':    
        
        w80min = 900
        w80max = 1450

        lvls = (1200,1300)
    if ID=='2QZJ':    
        
        w80min = 1600
        w80max = 2000
        
        lvls = (1900,1950)
        
    
    if ID=='LBQS':    
              
        w80min = 1300
        w80max = 2000
        
        lvls = (1800, 1900)
        
    
    cmap = plt.cm.coolwarm
    cmap.set_bad(color='black')
    
    fl = ax.imshow(smt, vmin=w80min, vmax=w80max, origin='low', cmap=cmap)
    
    ax.contour(smt, levels=lvls, colors='w', linestyles='dashed')
    
    axcbar0 = f.add_axes([ 0.58, 0.04 , 0.02, 0.23]) #plt.axes([0.055+ 0.2469,0.397,0.189,0.03])
    axcbar0.tick_params(direction='in')        
    cbar0 = f.colorbar(fl, cax=axcbar0 ,orientation='vertical')#, ticks= [min_flux,(min_flux+max_flux)/2, max_flux])
    #axcbar0.tick_params(axis='y',left='off',labelleft='off',right='off',labelright='off')
    #axcbar0.tick_params(axis='x',bottom='off',labelbottom='off',top='off',labeltop='off')
    
    axcbar0.set_ylabel('W80 (km/s)')
    
    ax.plot( cts[0], cts[1], 'r*', markersize=12)
    ax.set_xlim(coor[0]+lim[0],coor[0]+lim[1])
    ax.set_ylim(coor[1]+lim[2],coor[1]+lim[3])
   
    emplot.overide_axes_labels(f,ax,(lim_sc), labelx=1, labely=1, tickin=2)
    

    ax.plot( (cent[0]-step, cent[0]-step), (cent[1]+3*step, cent[1]-3*step), st, alpha=alp, linewidth=lnw)
    ax.plot( (cent[0]+step, cent[0]+step), (cent[1]+3*step, cent[1]-3*step), st, alpha=alp, linewidth=lnw)           
    ax.plot( (cent[0]-3*step, cent[0]-3*step), (cent[1]+3*step, cent[1]-3*step), st, alpha=alp, linewidth=lnw)
    ax.plot( (cent[0]+3*step, cent[0]+3*step), (cent[1]+3*step, cent[1]-3*step), st, alpha=alp, linewidth=lnw)
       
    ax.plot( (cent[0]-3*step, cent[0]+3*step), (cent[1]-3*step, cent[1]-3*step), st, alpha=alp, linewidth=lnw)
    ax.plot( (cent[0]-3*step, cent[0]+3*step), (cent[1]+1*step, cent[1]+1*step), st, alpha=alp, linewidth=lnw)     
    ax.plot( (cent[0]-3*step, cent[0]+3*step), (cent[1]-1*step, cent[1]-1*step), st, alpha=alp, linewidth=lnw)
    ax.plot( (cent[0]-3*step, cent[0]+3*step), (cent[1]+3*step, cent[1]+3*step), st, alpha=alp, linewidth=lnw)

# =============================================================================
#     Flux maps
# =============================================================================
    ax = f.add_axes([0.65, 0.04, 0.35, 0.23], projection=IFU_wcs)
    
    cmap = plt.cm.viridis
    cmap.set_bad(color='black')
    
   
    med = (Maps[4]*0.000188)*1e-13/(0.1*0.1)*1e14
    
    smt = smooth(med,1.5)
    use = np.isnan(med)
    smt[use] = np.nan
    
    fl = ax.imshow(med, origin='low', cmap=cmap)
    
    axcbar0 = f.add_axes([ 0.89, 0.04 , 0.02, 0.23]) #plt.axes([0.055+ 0.2469,0.397,0.189,0.03])
    axcbar0.tick_params(direction='in')        
    cbar0 = f.colorbar(fl, cax=axcbar0 ,orientation='vertical')
    
    axcbar0.set_ylabel('[OIII] SB \n($10^{-14}$ ergs/s/arcsec$^{2}$)')
    
    ax.set_xlim(coor[0]+lim[0],coor[0]+lim[1])
    ax.set_ylim(coor[1]+lim[2],coor[1]+lim[3])
    
    ax.plot( cts[0], cts[1], 'r*', markersize=12)
    
    emplot.overide_axes_labels(f,ax,(lim_sc), labelx=1, labely=1, tickin=2)
    
    ax.plot( (cent[0]-step, cent[0]-step), (cent[1]+3*step, cent[1]-3*step), st, alpha=alp, linewidth=lnw)
    ax.plot( (cent[0]+step, cent[0]+step), (cent[1]+3*step, cent[1]-3*step), st, alpha=alp, linewidth=lnw)           
    ax.plot( (cent[0]-3*step, cent[0]-3*step), (cent[1]+3*step, cent[1]-3*step), st, alpha=alp, linewidth=lnw)
    ax.plot( (cent[0]+3*step, cent[0]+3*step), (cent[1]+3*step, cent[1]-3*step), st, alpha=alp, linewidth=lnw)
       
    ax.plot( (cent[0]-3*step, cent[0]+3*step), (cent[1]-3*step, cent[1]-3*step), st, alpha=alp, linewidth=lnw)
    ax.plot( (cent[0]-3*step, cent[0]+3*step), (cent[1]+1*step, cent[1]+1*step), st, alpha=alp, linewidth=lnw)     
    ax.plot( (cent[0]-3*step, cent[0]+3*step), (cent[1]-1*step, cent[1]-1*step), st, alpha=alp, linewidth=lnw)
    ax.plot( (cent[0]-3*step, cent[0]+3*step), (cent[1]+3*step, cent[1]+3*step), st, alpha=alp, linewidth=lnw)

    emplot.new_ALMA_contour_plot(storage, ax, both=0, prj = '4QSO'   )
    
    f.savefig(PATH+'Four_Quasars/Graphs/Paper_plots/FigApp_'+ID+'_OIII_regions_wqso.pdf', bbox_inches = 'tight')
    
     


plt.show()
        
        
     
