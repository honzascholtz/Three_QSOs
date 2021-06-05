#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 30 15:04:31 2019

@author: jansen
"""

#importing modules
import numpy as np
import matplotlib.pyplot as plt ; plt.ioff()

from astropy.io import fits as pyfits
from astropy import wcs
from astropy.table import Table, join, vstack
from matplotlib.backends.backend_pdf import PdfPages
import pickle
from scipy.optimize import curve_fit

import Graph_setup as gst 
import glob
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
fsz = gst.graph_format(Labelsize=15)

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]

def find_nearest_idx(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

def gauss(x,k,sig):
    mean=0
    expo= -((x-mean)/(np.sqrt(2)*sig))**2    
    return k*e**expo

def straight(x, b):
    y = b*np.ones_like(x)
    return y 

import scipy.optimize as opt

from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.nddata import Cutout2D
import glob

from scipy.misc import imresize
import Plotting_tools as emplot

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]

def find_nearest_idx(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx



Sample = Table.read(PATH+'Four_Quasars/Four_Quasars.fits')
New = Sample




ylims = {'LBQS':( -2.4,2.4)}
ylims['HB89'] = (-2.5,3.2)
ylims['2QZJ'] = (-1,1.25)
ylims['XID_2028']= (-0.06,0.15)


lims_hal={'HB89': np.array([-9,14,-11,12])}
lims_hal['2QZJ'] = np.array([-8,10,-8,10])
lims_hal['LBQS'] = np.array([-10,10,-8,12])
lims_hal['XID_2028'] = np.array([-9,13,-9,13])


column = 0
row = 0

alph = 0.1
axes_width  = 0.20
axes_height = 0.2

im_width  = 0.15
im_height = 0.2

corner_x = 0.03
corner_y = 0.79
    

Fluxes = Table.read(PATH+'Four_Quasars/Catalogues/Target_table_total.fits')

fluxes = np.array([1.3, 1.5, 1.52,0.14])

h = plt.figure(figsize=(15,15))

axes_width = 0.165
corner_x = 0.05
corner_y = 0.88

row= 0 
column = 0


height = 0.1
height_off = 0.12
 
j=0

from scipy.interpolate import interp1d


Halpha_sizes = np.zeros(4)
Halpha_sizes_blr = np.zeros(4)

itere= np.array([2,1,0])


for i in itere:
    ID = Sample['ID'][i].strip()
    
# =============================================================================
#  IMAGES
# ==========================================================================

    if i==3:
        
        Ra_opt = Sample['RA'][i]
        Dec_opt = Sample['DEC'][i]
    
    else:
        Gaia = Table.read(PATH+'Four_Quasars/Catalogues/'+ID+'_Gaia.tbl' , format='ipac')
    
    
        Ra_opt = float(Gaia['ra'])
        Dec_opt = float(Gaia['dec'])
        
        
    new_size = 4.
   
    
    afile = glob.glob(PATH+'Four_Quasars/ALMA/Final_set/'+ID+'_l*')
    
    cont_file =glob.glob(PATH+'Four_Quasars/ALMA/Final_set/'+ID+'_l_clean.pbcor.*')[0]        
    cont= (pyfits.getdata(cont_file)[0,0,:,:])
       
    header = pyfits.getheader(cont_file)
    
    wcs_alma_header=  wcs.WCS(header)
    
    
    
    ###############################################################
    # Find the center of the image
    alm_world= np.array([[Ra_opt, Dec_opt,0,0]])    
    alm_pixcrd = wcs_alma_header.wcs_world2pix(alm_world, 0) # WCS transform
        
    alma_x= int(alm_pixcrd[0,0]) # X pixel
    alma_y= int(alm_pixcrd[0,1]) # Y pixel
        
    position = np.array([alma_x, alma_y])
                
    #print header
    pixscale = abs(header['CDELT1']*3600)
        
    cont_wcs= wcs.WCS(header).celestial
       
    cutout = Cutout2D(cont, position, new_size/pixscale, wcs=cont_wcs,mode='partial') 
        
    cont_c_d = cutout.data
    cont_c_w = cutout.wcs
    
    ax = h.add_axes((corner_x,corner_y-row*height_off, 0.10,height), projection=cont_c_w)
     
    rms = np.std(cont_c_d- cont_c_d.mean(axis=0))
    
    ax.imshow(cont_c_d,vmin=-rms, vmax=3*rms, origin='low')
    
    linew=2.
    ax.contour(cont_c_d,  transform= ax.get_transform(cont_c_w), colors='red', linestyles='solid' ,levels=( 2.0*rms ,3*rms,4*rms, 5*rms), alpha=0.5, linewidths=linew)
    
    ax.contour(cont_c_d,  transform= ax.get_transform(cont_c_w), colors='red', linestyles='dashed' ,levels=( -3.0*rms ,-2*rms,-1*rms, -0.5*rms), alpha=0.5, linewidths=linew)
    
    from matplotlib.patches import Circle
    
    lmn = new_size
    sz = header['BMAJ']/header['CDELT2']
    c = Circle((lmn*2, lmn*2), sz/2, edgecolor='firebrick', facecolor='firebrick', alpha=0.9)
     
    ax.add_patch(c)
    
    
    #ax.set_title(ID+ ' (%.2f' %(header['BMAJ']*3600)+'x %.2f' %(header['BMIN']*3600)+')')
    
    
    shapes = np.shape(cont_c_d)
    
    if ID=='HB89':
        IDp='HB8903'
    if ID=='2QZJ':
        IDp='2QZJ00'       
    if ID=='LBQS':
        IDp= 'LBQS01'
    
    ax.text(shapes[0]*0.05,0.8*shapes[0], IDp, color='w', fontsize=12)
    
    
    
    lim_sc = np.zeros(4)
    lim_sc[0] = - new_size/2
    lim_sc[1] = new_size/2
    lim_sc[2] = -new_size/2
    lim_sc[3] = new_size/2
    
    
    if (row==2):
        emplot.overide_axes_labels(h,ax,(lim_sc), labelx=1, labely=1,tickin=1, labelsize=15)       
    
    else:
        emplot.overide_axes_labels(h,ax,(lim_sc), labelx=0, labely=1,tickin=1, labelsize=15)
        
    
# =============================================================================
#     ALMA UV DIST
# =============================================================================
    print ('Starting ALMA uv')
    ax = h.add_axes((corner_x+0.157,corner_y-row*height_off, 0.17,height))
    
    ax.tick_params(direction='in')  
    
    Data = np.loadtxt('/Users/jansen/Google Drive/Astro/Four_Quasars/ALMA/UV_data/'+ID+'_phase_binnedvisibilities.txt', skiprows=0)    
    Data = Data[np.where((Data[:,0]>0)& (Data[:,0]<1500))[0],:]
    

    
    ax.tick_params(direction='in')  
    
    ax.errorbar(Data[:,0], Data[:,1], yerr=Data[:,2], fmt='o')
    #ax.set_title(ID)
    
    x = Data[:,0]
    y = Data[:,1]
    yr = Data[:,2]
    
    # Setting the            
    initial_guess = (1.0, 260.)
    
    popt, pcov = opt.curve_fit(gauss, x, y, p0=initial_guess, sigma=yr)
    er = np.sqrt(np.diag(pcov))
    xpl = np.linspace(min(x), 1500, 1000)
    ax.plot(xpl, gauss(xpl,*popt))
    
    
    chi2_gauss = (gauss(x,*popt)-y)**2*yr**-2
    chi2_gauss = sum(chi2_gauss)
    BIC_gauss = chi2_gauss + len(initial_guess)*len(y)
    

    print (ID, ' Flux ',  popt[0], ' +- ', er[0] )

    Size = 3600*(180./np.pi)*np.sqrt(2*np.log(2)/(np.pi*(popt[1]+er[1])*1000)**2  )
    size_erp = Size - 3600*(180./np.pi)*np.sqrt(2*np.log(2)/(np.pi*(popt[1]+er[1])*1000)**2  )
    size_erm =  3600*(180./np.pi)*np.sqrt(2*np.log(2)/(np.pi*(popt[1]-er[1])*1000)**2  )-Size
                        
    print (ID, ' Size ', Size, '+- ', size_erp, size_erm , ' arcseconds')
    
    
    initial_guess = ( 0.8)
    
    popt, pcov = opt.curve_fit(straight, x, y, p0=initial_guess, sigma=yr)
    er = np.sqrt(np.diag(pcov))
    xpl = np.linspace(min(x), 1500, 1000)
    ax.plot(xpl, straight(xpl,*popt), 'k--')
    
    print( ID, ' Flux point source ', popt[0])
    
    chi2_straight = (straight(x,*popt)-y)**2*yr**-2
    chi2_straight = sum(chi2_straight)
    BIC_straight = chi2_straight + 1*len(y)
    
    dBIC = BIC_straight-BIC_gauss
    print('Delta BIC, ', dBIC)
    print( ' ')
        
    #ax.set_ylim(lims[0], lims[1])
    ax.set_xlim(-5,950)
    lims = ylims[ID]
    ax.set_ylim(lims[0], lims[1])
    
    ax.text(300. , (lims[1]/2 ),r'r$_{e}$ = %.2f$\pm$'%(Size/2.354) +'%.2f arcsec' %(size_erm/2.354) )
    
    ax.set_ylabel('Amplitude (mJy)')
    
    if row==2:
        ax.set_xlabel(r'uv distance (k$\lambda$)')
    

    
    row+=1
    j+=1
    
#plt.savefig(PATH+'Four_Quasars/Graphs/Paper_plots/Fig2_ALMA_data_sizes.pdf', bbox_inches = 'tight')


#Halpha_sizes[0] = -Halpha_sizes_blr[0]
#Halpha_sizes[1] = -Halpha_sizes_blr[1]
#Halpha_sizes[2] = np.sqrt(Halpha_sizes_blr[2]**2 + Halpha_sizes_blr[2]**2 )


    
plt.show()