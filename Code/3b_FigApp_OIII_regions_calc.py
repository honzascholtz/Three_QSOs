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

import Graph_setup as gst 
import Tools_IFU as IFU
import Tools_plotting as emplot
import Tools_fitting as emfit
import Tools_path as ph
    
fsz = gst.graph_format()


binning = 'Nearest'



plot_it = 0
#0,4

Sample = Table.read(ph.MyPATH+'Four_Quasars.fits')



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


itere=np.array([2])
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
    
    storage_O = storage_H = IFU.load_cube(ph.MyPATH+'SINFONI/'+H_file+'.fits', z, ID, 'Sinfoni')

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
    
    ID = storage['X-ray ID']
    
    Save_spec = np.zeros((4,len(Spectra)))
    
    Save_spec[0,:] = wave
    Save_spec[1,:] = Spectra
    Save_spec[2,:] = error
    Save_spec[3,:] = storage['sky_clipped_1D']
    
    ID = storage['X-ray ID']
    
    np.savetxt(ph.MyPATH+'Results_storage/Spectrums/'+ID+'_OIII_inn.txt', Save_spec)
    
    hbw=4861
    if (ID=='HB89'):
        hbw=4870.
    if (ID=='LBQS'):
        hbw=4872
    outo = emfit.fitting_OIII_Hbeta_qso_mul(wave,Spectra, error,z, hbw=hbw)
    
    cmap = plt.cm.viridis
    cmap.set_bad(color='black')
    
    z = (outo.params['o3rn_center'].value/(5008./1e4))-1
    
    w80 = IFU.W80_mes(outo, 'OIII', 0)
    g, ax = plt.subplots(1)
    
    emplot.plotting_OIII_Hbeta(wave, Spectra, ax, outo, 'mul',z, title=0)
    
    Hb_flux = np.array([IFU.flux_measure_ind(outo, wave, 'Hb', use='BPT')])
    
    np.savetxt(ph.MyPATH+'Results_storage/'+ID+'_Hbeta_nuclear_flux.txt',Hb_flux)
    
    print('Nuclear Hb flux is ', Hb_flux)
    
    
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
#             Conitnuing on 
# =============================================================================
            flux = np.ma.array(data=storage['flux'].data, mask= mask) 
    
            Spectra = np.ma.sum(flux, axis=(1,2))
            Spectra = np.ma.array(data=Spectra.data, mask=storage['sky_clipped_1D'])
            error = IFU.STD_calc(wave/(1+z)*1e4, Spectra, 'OIII')* np.ones_like(Spectra)
            
            chi_list= np.array([])
            out_list= []
            
            #out = emfit.fitting_OIII_Hbeta_qso_mul(wave,Spectra, error,z, decompose=outo, offn=450, offw=700)
                        
            for ofn in np.array([0,50,100,150]):
                out = emfit.fitting_OIII_Hbeta_qso_mul(wave,Spectra, error,z, decompose=outo, offn=ofn)
                
                chi_list = np.append(chi_list, out.chisqr)
                out_list.append(out)
            
            out= out_list[np.argmin(chi_list)]
            print (np.array([-100,-50,0,50,100,150])[np.argmin(chi_list)])
            
            
            
            if (ID=='LBQS') & (iy==0):
                if ix==1:
                    out = emfit.fitting_OIII_Hbeta_qso_mul(wave,Spectra, error,z, decompose=outo, offn=50)
                if ix==2:
                    out = emfit.fitting_OIII_Hbeta_qso_mul(wave,Spectra, error,z, decompose=outo, offn=250, offw=700)
                else:
                    out = emfit.fitting_OIII_Hbeta_qso_mul(wave,Spectra, error,z, decompose=outo, offn=200, offw=700)
                    
            
                    
            #if ID=='HB89':
            #    out = emfit.fitting_OIII_Hbeta_qso_mul(wave,Spectra, error,z, decompose=outo, offn=-100)
            
            if ID=='2QZJ':
                if ix==2:
                    out = emfit.fitting_OIII_Hbeta_qso_mul(wave,Spectra, error,z, decompose=outo, offw=600, o3n=0)
                    
                    if iy==1:
                        out = emfit.fitting_OIII_Hbeta_qso_mul(wave,Spectra, error,z, decompose=outo, offw=700, o3n=1)
                    
                if (ix==0)&(iy==0):
                    out = emfit.fitting_OIII_Hbeta_qso_mul(wave,Spectra, error,z, decompose=outo, offw=600, o3n=0)
                    
                
            
            
            print('Whole region')
            print (out.params['o3rn_fwhm'].value/out.params['o3rn_center'].value*3e5)

            w80 = IFU.W80_mes(out, 'OIII', 0)
    
            emplot.plotting_OIII_Hbeta(wave, Spectra, ax, out, 'mul',z, title=0)
            print ('W80 ', w80[2])
            
            off = (outo.params['o3rn_center'].value - out.params['o3rw_center'].value)/outo.params['o3rn_center'].value*3e5
            print('Offset ', -off) 
            
            broadfwhm =out.params['o3rw_fwhm'].value/out.params['o3rw_center'].value*3e5
            print('broad fwhm ', broadfwhm) 
            ymin, ymax = ax.get_ylim()
            
            if ID=='2QZJ':
                if ix==2:
                    ymax = ymax*2
                    
                elif (ix==0)&(iy==0):
                    ymax = ymax*2
                
                else:
                    ymax = ymax*1.5
                    
                ax.set_ylim(ymin, ymax)
            
                    
            print (out.params['o3rn_height'].value)
            
            ax.axvline(outo.params['o3rn_center'].value/(1+z)*1e4, linestyle='--', color='k', alpha=0.3)
            
            F_broad = IFU.flux_measure_ind(out, wave, 'OIII',use='bro' )
                       
           
            
            ax.text(4720, ymax*0.85  , r'W80 %.0f km s$^{-1}$' %(w80[2]))
            ax.text(4720, ymax*0.75  , r'Offset %.0f km s$^{-1}$' %(-off))
            ax.text(4720, ymax*0.65  , r'Broad Flux %.1f $\times 10^{-16}$ ergs/s' %(F_broad*1e16))
            #ax.text(4720, ymax*0.1  , str(int(cent[0]+indy*step*2))+', '+ str(int(cent[1]+indx*step*2)) )
            #ax.text(4720, ymax*0.2  , str(indx)+', '+ str(indy) )
            
            
            ax.set_xlim(4700,5100)
            
            if indx==-1:
                ax.set_xlabel('Rest-frame wavelength (ang)', fontsize=12)
            
            
            if indy==-1:
                ax.set_ylabel(r'x$10^{-17}$ (erg/s/cm$^2$/ang)', fontsize=12)
            
            # Saving the data to use elsewhere 
            Save_spec = {'wave' : wave}
            Save_spec['data'] = Spectra
            Save_spec['error'] = error
            Save_spec['z'] = z
            Save_spec['total'] = out.eval(x=wave)
            
            Save_spec['cont'] = out.eval_components(x=wave)['linear']
            
            Save_spec['Hbw'] = out.eval_components(x=wave)['Hbw_']        
            Save_spec['Hbn'] = out.eval_components(x=wave)['Hbn_']
            
            Save_spec['o3rn'] = out.eval_components(x=wave)['o3rn_']
            Save_spec['o3bn'] = out.eval_components(x=wave)['o3bn_']
            
            Save_spec['o3rw'] = out.eval_components(x=wave)['o3rw_']
            Save_spec['o3bw'] = out.eval_components(x=wave)['o3bw_']
            
            Save_spec['W80'] = w80
            Save_spec['offset'] = off
            Save_spec['Broad flux'] = F_broad
            
                        
            #if (iy==0) & (ix==0): 
            #emplot.save_obj(Save_spec, 'Results_storage/Spectrums/Regions/'+ID+'_OIII_'+str(iy)+'_'+str(ix))


plt.show()
        
        
     
