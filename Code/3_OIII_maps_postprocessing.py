#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  3 10:19:30 2020

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



center = {'LBQS' : np.array([53,47])}
center['HB89'] = np.array([53.7,47.9])
center['2QZJ'] = np.array([53.0, 47.5])

lims={'HB89':  np.array([-10,10,-10,10])}
lims['2QZJ'] = np.array([-10,12,-10,12])
lims['LBQS'] = np.array([-10,10,-10,10])


def smooth(image,sm):
    
    from astropy.convolution import Gaussian2DKernel
    
    from astropy.convolution import convolve
    
    gauss_kernel = Gaussian2DKernel(sm)

    con_im = convolve(image, gauss_kernel)
    
    #con_im = con_im#*image/image
    
    return con_im  
def post_process(Results):
    # VLT SINFONI Wavelength setup
    wave_obs = np.linspace(1.43549998, 1.85, 2202)
    
    # It is a way to create the model and then "fit it" to create the right results class
    from lmfit.models import GaussianModel
    gmodel = GaussianModel(prefix='o3rn_') + GaussianModel(prefix='o3rw_')
    
    
    params = gmodel.make_params(o3rn_center=Results[0], o3rn_sigma =Results[1], o3rn_amplitude=Results[2], \
                                o3rw_center=Results[3], o3rw_sigma =Results[4], o3rw_amplitude=Results[5])

    y_eval = gmodel.eval(params, x=wave_obs)
    out = gmodel.fit(y_eval, params, x=wave_obs)
    
    # Calculates v10, v90, W80 and v50
    w80 = IFU.W80_mes(out, 'OIII', 0)
        
    peak = wave_obs[np.argmax(y_eval)]
    flux = sum(y_eval)
        
    return w80[0], w80[1], w80[2] , peak, flux, w80[3], y_eval



Sample = Table.read(ph.MyPATH+'Four_Quasars.fits')


itere=np.array([2])
for i in itere:

    ID = Sample['ID'][i].strip()
    z= Sample['z'][i]

    print(ID, z)
    
    Results = pyfits.getdata(ph.MyPATH+'Results_storage/OIII/'+ID+'_Individual_spaxel_fit.fits')
    prhdr = pyfits.getheader(ph.MyPATH+'Results_storage/OIII/'+ID+'_Individual_spaxel_fit.fits')   
    
    
    
    Maps = True
    if Maps==True:
        shapes = np.shape(Results)
        
        Maps = np.zeros((6,shapes[1], shapes[2]))
        Maps[:,:,:] = np.nan
        
        wave_obs = np.linspace(1.43549998, 1.85, 2202)
        
        Cube = np.zeros((2202,shapes[1], shapes[2]))
        
        
        for i in range(shapes[1]):
            if Results[0,i,50] >0:
                print (i)
            for j in range(shapes[2]):
                if Results[0,i,j] >0:
                    res = post_process(Results[:,i,j])
                    
                    l = res[:6]
                    
                    Maps[:,i,j] = l                     
                    Cube[:,i,j] = res[6]
        
        
        hdu = pyfits.PrimaryHDU(Maps, header=prhdr)
        hdulist = pyfits.HDUList([hdu])  
        hdulist.writeto(ph.MyPATH+'Results_storage/'+ID+'_OIII_mapsind.fits', overwrite=True)
        
        hdu = pyfits.PrimaryHDU(Cube, header=prhdr)
        hdulist = pyfits.HDUList([hdu])  
        hdulist.writeto(ph.MyPATH+'Results_storage/'+ID+'_OIII_modelind.fits', overwrite=True)
    
    Maps = pyfits.getdata(ph.MyPATH+'Results_storage/'+ID+'_OIII_mapsind.fits')
    # =============================================================================
    # Plotting stuff
    # =============================================================================
    if ID=='HB89':    
        v90min = 200
        v90max = 550
        
        v10min = -700
        v10max = -250
        
        w80min = 1000
        w80max = 1600
        
        v50min = -200
        v50max = 50
        
    if ID=='2QZJ':    
        v90min = 100
        v90max = 400
        
        v10min = -2000
        v10max = -1000
        
        w80min = 1200#/2.354
        w80max = 2200#/2.354
        
        v50min = 10
        v50max = -700
        #Maps[2,:,:] = Maps[2,:,:]/2.354
    
    if ID=='LBQS':    
        v90min = 600
        v90max = 900
        
        v10min = -1100
        v10max = -300
        
        w80min = 1100
        w80max = 2000
        
        v50min = 200
        v50max = -100
        
        
    
    lim = lims[ID]
    
    f = plt.figure(figsize=(15,9))
    # =============================================================================
    # OIII Flux
    ax1 = f.add_axes([0.05, 0.55, 0.3, 0.4])
    ax1.set_title('OIII flux')
    ax1.imshow(Maps[4,:,:], origin='low')
    
    
    
    cts=coor = center[ID]
    ax1.plot( cts[0], cts[1], 'r*', markersize=12)
    
    ax1.set_xlim(coor[0]+lim[0],coor[0]+lim[1])
    ax1.set_ylim(coor[1]+lim[2],coor[1]+lim[3])
    
    # OIII W80
    ax = f.add_axes([0.37, 0.55, 0.3, 0.4])
    
    if ID=='2QZJ':
        ax.set_title(r'$\sigma$ vel (2.354 conversion)')
        ax.set_title('w80 vel')
    else:
        ax.set_title('w80 vel')
    fl = ax.imshow(Maps[2,:,:], vmin=w80min, vmax=w80max, origin='low')
    
    axcbar0 = plt.axes([ 0.63, 0.55 , 0.02, 0.4]) #plt.axes([0.055+ 0.2469,0.397,0.189,0.03])
    axcbar0.tick_params(direction='in')        
    cbar0 = f.colorbar(fl, cax=axcbar0 ,orientation='vertical')#, ticks= [min_flux,(min_flux+max_flux)/2, max_flux])
    #axcbar0.tick_params(axis='y',left='off',labelleft='off',right='off',labelright='off')
    #axcbar0.tick_params(axis='x',bottom='off',labelbottom='off',top='off',labeltop='off')
    
    cts = center[ID]
    ax.plot( cts[0], cts[1], 'r*', markersize=12)
    
    if ID=='2QZJ':
        out = Maps[2,:,:].copy()
        out[np.isnan(out)] = 0
        
        lvls = np.array([1800,2000])
        
        
        hdu = pyfits.PrimaryHDU(out, header=prhdr)
        hdulist = pyfits.HDUList([hdu])  
        hdulist.writeto(ph.MyPATH+'Results_storage/2QZJ_outflow.fits', overwrite=True)
        
        np.savetxt(ph.MyPATH+'Results_storage/2QZJ_outflow_cont.txt',lvls)
        
        
        ax1.contour(out, levels=lvls,linestyles='dashed', colors='red')
        
    
    if ID=='HB89':
        out = Maps[2,:,:].copy()
        out[np.isnan(out)] = 0
        out[55:,:] =0
        out[:40,:] =0
        
        lvls= (1000,1200)
        
        
        ax.contour(out, levels=lvls,linestyles='dashed', colors='red')
    
    ax.set_xlim(coor[0]+lim[0],coor[0]+lim[1])
    ax.set_ylim(coor[1]+lim[2],coor[1]+lim[3])
    
    # =============================================================================
    # BLR Hbeta 
    ax = f.add_axes([0.69, 0.55, 0.3, 0.4])
    ax.set_title('BLR Hbeta')
    ax.imshow(Results[8,:,:], origin='low')
    cts = center[ID]
    ax.plot( cts[0], cts[1], 'r*', markersize=12)
    
    coor=cts
    
    ax.set_xlim(coor[0]+lim[0],coor[0]+lim[1])
    ax.set_ylim(coor[1]+lim[2],coor[1]+lim[3])
    
    # =============================================================================
    # v10 Map
    ax = f.add_axes([0.05, 0.05, 0.3, 0.4])
    ax.set_title('v10 vel')
    fl = ax.imshow(Maps[0,:,:],vmin=v10min, vmax=v10max, origin='low')
    
    axcbar0 = plt.axes([ 0.3, 0.05 , 0.02, 0.4]) #plt.axes([0.055+ 0.2469,0.397,0.189,0.03])
    axcbar0.tick_params(direction='in')        
    cbar0 = f.colorbar(fl, cax=axcbar0 ,orientation='vertical')
    
    ax.set_xlim(coor[0]+lim[0],coor[0]+lim[1])
    ax.set_ylim(coor[1]+lim[2],coor[1]+lim[3])
    
    cts = center[ID]
    ax.plot( cts[0], cts[1], 'r*', markersize=12)
    
    if ID=='HB89':
        out = Maps[0,:,:].copy()
        #out[np.isnan(out)] = 0
        out[55:,:] =0
        out[:40,:] =0
        
        lvls = np.array([-700,-550])
        
        hdu = pyfits.PrimaryHDU(out, header=prhdr)
        hdulist = pyfits.HDUList([hdu])  
        hdulist.writeto(ph.MyPATH+'Results_storage/HB89_outflow.fits', overwrite=True)
        
        np.savetxt(ph.MyPATH+'Results_storage/HB89_outflow_cont.txt',lvls)
        
        ax1.contour(out, levels=lvls,linestyles='dashed', colors='red')
    
    # =============================================================================
    # v90 Map 
    ax = f.add_axes([0.37, 0.05, 0.3, 0.4])
    ax.set_title('v90 vel')
    fl = ax.imshow(Maps[1,:,:],vmin=v90min, vmax=v90max, origin='low')
    
    axcbar0 = plt.axes([ 0.63, 0.05 , 0.02, 0.4]) #plt.axes([0.055+ 0.2469,0.397,0.189,0.03])
    axcbar0.tick_params(direction='in')        
    cbar0 = f.colorbar(fl, cax=axcbar0 ,orientation='vertical')
    
    ax.set_xlim(coor[0]+lim[0],coor[0]+lim[1])
    ax.set_ylim(coor[1]+lim[2],coor[1]+lim[3])
    
    cts = center[ID]
    ax.plot( cts[0], cts[1], 'r*', markersize=12)
    
    # =============================================================================
    # v50 Map 
    ax = f.add_axes([0.69, 0.05, 0.3, 0.4])
    ax.set_title('v50')
    fl = ax.imshow(Maps[5,:,:], vmin=v50min, vmax=v50max , origin='low')
    
    axcbar0 = plt.axes([ 0.94, 0.05 , 0.02, 0.4]) #plt.axes([0.055+ 0.2469,0.397,0.189,0.03])
    axcbar0.tick_params(direction='in')        
    cbar0 = f.colorbar(fl, cax=axcbar0 ,orientation='vertical')
    
    ax.set_xlim(coor[0]+lim[0],coor[0]+lim[1])
    ax.set_ylim(coor[1]+lim[2],coor[1]+lim[3])
    
    cts = center[ID]
    ax.plot( cts[0], cts[1], 'r*', markersize=12)
    
    f, ax = plt.subplots(1)
    
    
    cnt =  5008.*(1+z)/1e4
    med = (Results[3]-cnt)/cnt*3e5
    
    if ID=='LBQS':
        vn = -200
        vx = 300
        
    if ID=='HB89':
        vn = -500
        vx = 100
    
    if ID=='2QZJ':
        vn = -300
        vx = -200
   
    fl  =ax.imshow(med,vmin=vn, vmax=vx, origin='low')
    
    plt.colorbar(fl)
    #ax.contour(med, levels=(-300,-275,-250), colors='red')
    
    ax.set_title('Broad offset')
    ax.set_xlim(coor[0]+lim[0],coor[0]+lim[1])
    ax.set_ylim(coor[1]+lim[2],coor[1]+lim[3])
    
    #plt.savefig(PATH+'Four_Quasars/Graphs/'+ID+'_OIII_maps.pdf')

plt.show()