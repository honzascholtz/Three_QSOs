#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 14 15:22:53 2020

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

def twoD_Gaussian(dm, amplitude, xo, yo, sigma_x, sigma_y, theta, offset): 
    x, y = dm
    xo = float(xo)
    yo = float(yo)    
    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
    g = offset + amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo) 
                            + c*((y-yo)**2)))
    return g.ravel()

def smooth(image,sm):
    
    from astropy.convolution import Gaussian2DKernel
    from scipy.signal import convolve as scipy_convolve
    from astropy.convolution import convolve
    
    gauss_kernel = Gaussian2DKernel(sm)

    con_im = convolve(image, gauss_kernel)
    
    #con_im = con_im#*image/image
    
    return con_im  

import Plotting_tools as emplot

lims={'HB89': np.array([-9,14,-11,12])}
lims['2QZJ'] = np.array([-20,20,-20,20])
lims['LBQS'] = np.array([-10,10,-8,12])



ID='2QZJ'


Image_narrow = pyfits.getdata(PATH+'KMOS_SIN/Results_storage/Halpha/'+ID+'_Nearest_spaxel_fit.fits')
IFU_header = pyfits.getheader(PATH+'KMOS_SIN/Results_storage/Halpha/'+ID+'_Nearest_spaxel_fit.fits')   
     
IFU_wcs= wcs.WCS(IFU_header).celestial
     
    

nrmm = Image_narrow[5].copy()
nrmm[np.isnan(nrmm)] = 0
nrm = nrmm.max()


nrxm  = Image_narrow[0].copy()
nrxm[np.isnan(nrxm)] = 0
nrx = nrxm.max()

shapes = np.shape(Image_narrow)

x = np.linspace(0, shapes[2]-1, shapes[2])
y = np.linspace(0, shapes[1]-1, shapes[1])
x, y = np.meshgrid(x, y)   

dm = (x,y)
poptb, pcov = curve_fit(twoD_Gaussian, dm, nrmm.ravel(), p0=(nrm, shapes[2]/2,shapes[2]/2, 5,5,0,0))

poptn, pcov = curve_fit(twoD_Gaussian, dm, nrxm.ravel(), p0=(nrx, shapes[2]/2,shapes[2]/2, 5,5,0,0))

plt.figure()

plt.imshow(nrmm,origin='low')

plt.plot(poptb[1], poptb[2], 'ro')


# =============================================================================
# Plotting
# =============================================================================
f = plt.figure()

ax= f.add_axes((0.1,0.1, 0.85,0.85), projection= IFU_wcs)

image = Image_narrow[0]- twoD_Gaussian(dm, poptn[0], poptn[1], poptn[2], poptb[3], poptb[4],poptb[5], poptn[6]  ).reshape(shapes[1], shapes[2])
image = smooth(image,1.)

use = np.isnan(Image_narrow[0])
image[use] = np.nan

ax.plot(poptb[1], poptb[2], 'bo')
ax.plot(poptn[1], poptn[2], 'ro')
ax.imshow( image , origin='low', vmin=0)

lim = lims[ID]

coor = np.array([poptb[1], poptb[2]])
ax.set_xlim(coor[0]+lim[0],coor[0]+lim[1])
ax.set_ylim(coor[1]+lim[2],coor[1]+lim[3])
    
    
deg_per_pix = IFU_header['CDELT2']
arc_per_pix = deg_per_pix*3600
lim_sc = lim*arc_per_pix
    
#emplot.overide_axes_labels(f,ax,(lim_sc), labelx=1, labely=1,tickin=1)

np.savetxt(PATH+'Four_Quasars/Results_storage/'+ID+'_Halpha_PSF_sub_smt.txt', image)


plt.show()