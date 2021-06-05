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
    wave_obs = np.linspace(1.43549998, 1.85, 2202)
    
    from lmfit.models import GaussianModel
    gmodel = GaussianModel(prefix='o3rn_') + GaussianModel(prefix='o3rw_')
    
    
    params = gmodel.make_params(o3rn_center=Results[0], o3rn_sigma =Results[1], o3rn_amplitude=Results[2], \
                                o3rw_center=Results[3], o3rw_sigma =Results[4], o3rw_amplitude=Results[5])

    y_eval = gmodel.eval(params, x=wave_obs)
    out = gmodel.fit(y_eval, params, x=wave_obs)
    
    
  
    
    
    fl_nar =  IFU.flux_measure_ind(out,wave_obs, 'OIII', use='nar', error=0)
    fl_bro =  IFU.flux_measure_ind(out,wave_obs, 'OIII', use='bro', error=0)
    
    return fl_nar, fl_bro        
    



Sample = Table.read(ph.MyPATH+'Four_Quasars.fits')


itere=np.array([0,1,2])
for i in itere:

    ID = Sample['ID'][i].strip()
    z= Sample['z'][i]

    print(ID, z)
    
    Results = pyfits.getdata(ph.MyPATH+'Results_storage/OIII/'+ID+'_Individual_spaxel_fit.fits')
    prhdr = pyfits.getheader(ph.MyPATH+'Results_storage/OIII/'+ID+'_Individual_spaxel_fit.fits')   
    
    
    
    Maps = True
    if Maps==True:
        shapes = np.shape(Results)
        
        Maps = np.zeros((2,shapes[1], shapes[2]))
        Maps[:,:,:] = np.nan
        
        wave_obs = np.linspace(1.43549998, 1.85, 2202)
        
        
        
        for i in range(shapes[1]):
            if Results[0,i,50] >0:
                print (i)
            for j in range(shapes[2]):
                if Results[0,i,j] >0:
                    res = post_process(Results[:,i,j])
                    
                    l = res
                    
                    Maps[:,i,j] = l                     
                    
        
        
        hdu = pyfits.PrimaryHDU(Maps, header=prhdr)
        hdulist = pyfits.HDUList([hdu])  
        hdulist.writeto(ph.MyPATH+'Results_storage/'+ID+'_OIII_mapsind_flux.fits', overwrite=True)
        
            
    Maps = pyfits.getdata(ph.MyPATH+'Results_storage/'+ID+'_OIII_mapsind_flux.fits')
    # =============================================================================
    # Plotting stuff
    # =============================================================================
    
    
    f, (ax1,ax2) = plt.subplots(1,2)
    
    ax1.imshow(Maps[0], origin='low')
    
    ax2.imshow(Maps[1], origin='low')
    

import Plotting_tools as emplot


f= plt.figure( figsize=(7,10))    
    
itere=np.array([2,1,0])
row=0
for i in itere:
    ID = Sample['ID'][i].strip()
    
    Maps = pyfits.getdata(ph.MyPATH+'Results_storage/'+ID+'_OIII_mapsind_flux.fits')
    Header = pyfits.getheader(ph.MyPATH+'Results_storage/'+ID+'_OIII_mapsind_flux.fits')
   
    
    IFU_wcs = wcs.WCS(Header).celestial
    
    ax= f.add_axes((0.13,0.70-(row*0.32), 0.35,0.28), projection= IFU_wcs)
    
    lim = lims[ID]
    coor = center[ID]
    
    ax.imshow(Maps[0], origin='low')
    
    ax.set_title(ID+' Nar')
                  
    ax.plot( coor[0], coor[1], 'r*', markersize=12)               
    
    ax.set_xlim(coor[0]+lim[0],coor[0]+lim[1])
    ax.set_ylim(coor[1]+lim[2],coor[1]+lim[3])
    
    deg_per_pix = Header['CDELT2']
    arc_per_pix = deg_per_pix*3600
    
    lim_sc = lim*arc_per_pix
    
    emplot.overide_axes_labels(f,ax,(lim_sc), labelx=1, labely=1,tickin=1)
    
    
    ax= f.add_axes((0.6,0.70-(row*0.32), 0.35,0.28), projection= IFU_wcs)
    
    lim = lims[ID]
    coor = center[ID]
    
    ax.imshow(Maps[1], origin='low')
                  
    ax.set_title(ID+' Bro')              
    
    ax.set_xlim(coor[0]+lim[0],coor[0]+lim[1])
    ax.set_ylim(coor[1]+lim[2],coor[1]+lim[3])
    
    deg_per_pix = Header['CDELT2']
    arc_per_pix = deg_per_pix*3600
    
    ax.plot( coor[0], coor[1], 'r*', markersize=12) 
            
    lim_sc = lim*arc_per_pix
    
    emplot.overide_axes_labels(f,ax,(lim_sc), labelx=1, labely=1,tickin=1)
    
    
    row+=1
    
    
    
    
plt.savefig(ph.MyPATH+'Graphs/OIII_flux.pdf', bbox_inches='tight')
    
for i in itere:
    ID = Sample['ID'][i].strip()
    
    Maps = pyfits.getdata(ph.MyPATH+'Results_storage/'+ID+'_OIII_mapsind_flux.fits')
    Header = pyfits.getheader(ph.MyPATH+'Results_storage/'+ID+'_OIII_mapsind_flux.fits')
    plt.figure()
    
    plt.imshow( Maps[0]/Maps[1],vmax=2, origin='low')      
    
    
plt.show()