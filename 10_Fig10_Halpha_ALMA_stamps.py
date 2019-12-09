#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 18 12:04:28 2018

@author: jansen
"""

#importing modules
import numpy as np
import matplotlib.pyplot as plt; plt.ioff()

from astropy.io import fits as pyfits
from astropy import wcs
from astropy.table import Table, join, vstack
from astropy.coordinates import SkyCoord

from matplotlib.backends.backend_pdf import PdfPages
import pickle
from scipy.optimize import curve_fit
import glob


import Graph_setup as gst 
import IFU_tools as IFU
import Plotting_tools as emplot

from scipy.misc import imresize


def smooth(image,sm):
    
    from astropy.convolution import Gaussian2DKernel
    from scipy.signal import convolve as scipy_convolve
    from astropy.convolution import convolve
    
    gauss_kernel = Gaussian2DKernel(sm)

    con_im = convolve(image, gauss_kernel)
    
    #con_im = con_im#*image/image
    
    return con_im  
 
    


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


New= Table.read(PATH+'Four_Quasars/Four_Quasars.fits')

    

def twoD_Gaussian((x, y), amplitude, xo, yo, sigma_x, sigma_y, theta, offset):
    xo = float(xo)
    yo = float(yo)    
    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
    g = offset + amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo) 
                            + c*((y-yo)**2)))
    return g.ravel()



###############################################################################
# Halpha nad ALMA
###############################################################################
k=1

column = 0
row = 0

lims={'HB89': np.array([-9,14,-11,12])}
lims['2QZJ'] = np.array([-20,20,-20,20])
lims['LBQS'] = np.array([-10,10,-8,12])


lims['XID_2028'] = np.array([-9,13,-9,13])



IFU_psfs = {'HB89': np.array([-1.1,-0.6,-0.08])/3600}
IFU_psfs['2QZJ'] = np.array([-0.8,-0.7, -0.1 ])/3600
IFU_psfs['LBQS'] = np.array([ 0.0,-0.55, -0.1])/3600

IFU_psfs['lid_1565'] = np.array([-1.3, 0.7, 0])/3600

IFU_psfs['XID_614'] = np.array([ 0.2,-0.9,-0.1])/3600
IFU_psfs['XID_751'] = np.array([-0.8,-1.6, 0.1])/3600
IFU_psfs['ALESS_75'] = np.array([-0.8,0.8, +0.04])/3600



# =============================================================================
# Halpha, ALMA and and Outflow map
# =============================================================================
#f = plt.figure(figsize=(8.3,8.3))
f = plt.figure(figsize=(6,14))

itere = np.array([2, 1,0])

for i in itere:
    ID = New['ID'][i].strip()
      
    z = New['z'][i]    
    
    print i,ID, z
    # Halpha    
    

    if ID=='XID_2028':
        
        Image_narrow = pyfits.getdata(PATH+'KMOS_SIN/Results_storage/Halpha/lid_1565_Nearest_spaxel_fit_sin.fits')
        IFU_header = pyfits.getheader(PATH+'KMOS_SIN/Results_storage/Halpha/lid_1565_Nearest_spaxel_fit.fits')   
        
        IFU_header['CRPIX1'] = 48.3
        IFU_header['CRPIX2'] = 46.4
        
        IFU_header['CRVAL1'] = 150.5470381497013
        IFU_header['CRVAL2'] = 1.618527398264767
    
    else:
        Image_narrow = pyfits.getdata(PATH+'KMOS_SIN/Results_storage/Halpha/'+ID+'_Nearest_spaxel_fit.fits')
        IFU_header = pyfits.getheader(PATH+'KMOS_SIN/Results_storage/Halpha/'+ID+'_Nearest_spaxel_fit.fits')   
        
    Hal_map = Image_narrow[0]
    
    if ID=='XID_2028':
        
        storage = emplot.load_obj('KMOS_SIN/Results_storage/Props/lid_1565_H')
    else:
        storage = emplot.load_obj('KMOS_SIN/Results_storage/Props/'+ID+'_H')
    
    storage['Cat_entry'] = New[i]
    storage['X-ray ID'] = ID
    
    shapes = np.shape(Hal_map)
    
    data = smooth(Hal_map,1)
    
    
    
    data[np.isnan(data)] = 0   
    
    if ID=='XID_2028':
        data[:42, :] =  0
        data[55:, :] =  0
        data[:, 53:] =  0
    
    loc = np.ma.where(data > 0.98*data.max())   
    

    p_x = np.median(loc[1])
    p_y = np.median(loc[0])
    
   
    
        
    IFU_wcs= wcs.WCS(IFU_header).celestial
    
    
    #ax= f.add_axes((0.09+(column*0.48),0.52-(row*0.46), 0.42,0.42), projection= IFU_wcs)
    ax= f.add_axes((0.15,0.69-(row*0.33), 0.7,0.3), projection= IFU_wcs)
    
    cmap = plt.cm.viridis
    cmap.set_bad(color='black')
    
    if ID=='XID_2028':
        ax.imshow(Hal_map, cmap=cmap, origin='low',vmax=5e-17, aspect='auto', interpolation='nearest')
        
    elif ID=='2QZJ':
        ax.imshow(Hal_map, cmap=cmap, origin='low',vmax=7e-18, aspect='auto', interpolation='nearest')
    
    else:
        ax.imshow(Hal_map, cmap=cmap, origin='low', aspect='auto', interpolation='nearest')
    
    coor = storage['Cube_cent'][1:3].copy()
    coor = np.array(coor, dtype=int)
    lim = lims[ID]
        
    ax.set_xlim(coor[0]+lim[0],coor[0]+lim[1])
    ax.set_ylim(coor[1]+lim[2],coor[1]+lim[3])
        
    ax.set_autoscale_on(False)
    
    ax.scatter(p_x,p_y,marker='o', color='red')
    
    arc = np.round(1./(IFU_header['CDELT2']*3600))
    step = int(arc/4)   
    cent = np.round(storage['Cube_cent'][1:3]).copy()
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
    
    
    
    emplot.new_ALMA_contour_plot(storage, ax, both=0, prj = '4QSO'   )
    
    ax.plot(coor[0], coor[1], 'r*', markersize=10)
    
    
    
 # =============================================================================
#     ALMA Position
# =============================================================================
    Img = pyfits.getdata(PATH+'Four_Quasars/ALMA/Final_set/'+ID+'_l_clean.pbcor.fits')[0,0,:,:]
    ALM_header = pyfits.getheader(PATH+'Four_Quasars/ALMA/Final_set/'+ID+'_l_clean.pbcor.fits')  
    
    Ra= New['RA'][i]
    Dec = New['DEC'][i]
    
    ALM_wcs= wcs.WCS(ALM_header).celestial
    
    # Finding the position of the Galaxy in pix scale
    opt_world= np.array([[Ra,Dec]])
    opt_pixcrd = ALM_wcs.wcs_world2pix(opt_world, 0) # WCS transform
    opt_x= (opt_pixcrd[0,0]) # X pixel
    opt_y= (opt_pixcrd[0,1]) # Y pixel
    
    Img[:int(opt_x-30),:] = 0
    Img[int(opt_x+30):,:] = 0
    
    Img[:,:int(opt_x-30)] = 0
    Img[:,int(opt_x+30):] = 0
    
    loc = np.ma.where(Img > 0.98*Img.max())
    
    p_x = np.median(loc[1])
    p_y = np.median(loc[0])
    
    ALM_wcs= wcs.WCS(ALM_header).celestial
    
    alm_world = ALM_wcs.wcs_pix2world(np.array([  [p_x, p_y  ] ]),0)
    alm_world = alm_world[0,:]
    
    ax.scatter(alm_world[0],alm_world[1],transform=ax.get_transform('fk5'), marker='o', color='blue')
    
    
    
    
    dx = lim[1] - lim[0]
    dy = lim[3] - lim[2]
    x = coor[0]+lim[0]
    y = coor[1]+lim[2]
    
      
# =============================================================================
#     if ID=='ALESS_75':
#         ax.text(dx*0.05+x, dy*0.9+y, ID+'.1', color='white')
#     
#     elif ID=='lid_1565':
#         ax.text(dx*0.05+x, dy*0.9+y, 'XID_2028', color='white')
#     
#     
#     else:
#         ax.text(dx*0.05+x, dy*0.9+y, ID, color='white')
# =============================================================================
    ax.text(dx*0.05+x, dy*0.9+y, ID, color='white')
    
    deg_per_pix = IFU_header['CDELT2']
    arc_per_pix = deg_per_pix*3600
    
    lim_sc = lim*arc_per_pix
    
    if row==2:
        emplot.overide_axes_labels(f,ax,(lim_sc), labelx=1, labely=1,tickin=1)

    else:
        emplot.overide_axes_labels(f,ax,(lim_sc), labelx=0, labely=1,tickin=1)
    row+=1
    '''
    if column==0:
        emplot.overide_axes_labels(f,ax,(lim_sc), labelx=0, labely=1,tickin=1)
    
    else:
        emplot.overide_axes_labels(f,ax,(lim_sc), labelx=0, labely=0,tickin=1)
    
    if row==2:
        emplot.overide_axes_labels(f,ax,(lim_sc), labelx=1, labely=0,tickin=1)
    
    else:
        emplot.overide_axes_labels(f,ax,(lim_sc), labelx=0, labely=0,tickin=1)
    
    if (row==1) & (column==0):
        emplot.overide_axes_labels(f,ax,(lim_sc), labelx=1, labely=1,tickin=1)
    
    if (row==1) & (column==1):
        emplot.overide_axes_labels(f,ax,(lim_sc), labelx=1, labely=0,tickin=1)
    
    
    
    
    if column==2:
        row+=1
        column=0
    '''
    k+=1






'''
ax.plot((-100),(-100), 'r--', label= 'FIR continuum \n- IFU matched')
ax.plot((-100),(-100), 'r-', label= 'FIR continuum  \n- High res.')

ax.plot((-100),(-100), 'bo', label= r'Centre of H$\alpha$')
ax.plot((-100),(-100), 'ro', label= r'Centre of FIR')

legend = f.legend(bbox_to_anchor=(0.73, 0.3), loc=2, borderaxespad=0.,
                  title=r'Image - Narrow H$\alpha$', fontsize='12' ,framealpha=1., facecolor='black')

plt.setp(legend.get_title(),fontsize='12')

for text in legend.get_texts():
    plt.setp(text, color = 'w')

legend.get_title().set_color('w')

'''  
f.savefig(PATH+'Four_Quasars/Graphs/Paper_plots/Halpha_ALMA.pdf', bbox_inches = 'tight')
plt.show()
