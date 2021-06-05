#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  9 22:41:55 2019

@author: jansen
"""

#importing modules
import numpy as np
import matplotlib.pyplot as plt

from astropy.io import fits as pyfits
from astropy import wcs
from astropy.table import Table, join, vstack
from matplotlib.backends.backend_pdf import PdfPages
import pickle
from scipy.optimize import curve_fit

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

import Graph_setup as gst 
import Tools_IFU as IFU
import Tools_plotting as emplot
import Tools_fitting as emfit
import Tools_path as ph

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


from astropy.modeling.models import Sersic2D
import matplotlib.pyplot as plt

x,y = np.meshgrid(np.arange(100), np.arange(100))

mod = Sersic2D(amplitude = 1, r_eff = 15, n=1, x_0=50, y_0=50,
               ellip=0.01, theta=0)



x,y = np.meshgrid(np.arange(100), np.arange(100))


noise = np.random.normal(0, 0.05, 10000).reshape(100,100)

galaxy = twoD_Gaussian((x,y), 1., 50,50, 11,11,0,0).reshape(100,100)+noise
galaxy_hole = galaxy.copy()
galaxy_bar = galaxy.copy()

x = np.linspace(50-9,50+9, 19)
x = np.array(x, dtype=int)

y = np.linspace(50-9,50+9, 19)
y = np.array(y, dtype=int)

for ix in x:
    for iy in y:
        
        r = np.sqrt((ix-50)**2+(iy-50)**2)
        
        if r<5:
            galaxy_hole[iy,ix] = galaxy[ix,iy]  - twoD_Gaussian((ix,iy), 1., 50,50, 11,11,0,0)#.reshape(100,100)


y = np.linspace(50-8,50+8, 17)
y = np.array(y, dtype=int)

x = np.linspace(47,99, 53)
x = np.array(x, dtype=int)

for ix in x:
    for iy in y:
        
        galaxy_bar[iy,ix] = galaxy[ix,iy]  - twoD_Gaussian((ix,iy), 1., 50,50, 11,11,0,0)#.reshape(100,100)




plt.imshow(galaxy_bar, origin='low')




# =============================================================================
# COG
# =============================================================================


rads = np.linspace(3,27,25)
fl_mod = np.zeros(25)
fl_mod_gap = np.zeros(25)
fl_mod_wei = np.zeros(25)



mask = np.ma.masked_invalid(galaxy).mask
mask[:,:] = True



for i in range(25):
    mask[:,:] = True
    for ix in range(100):
        for iy in range(100):
            dist = np.sqrt((ix- 50)**2+ (iy- 50)**2)
            if dist< i:
                mask[ix,iy] = False
    
    
    fl_mod[i] = np.ma.sum(np.ma.array(data=galaxy, mask=mask))
    
    fl_mod_gap[i] = np.ma.sum(np.ma.array(data=galaxy_hole, mask=mask))
    
    fl_mod_wei[i] = np.ma.sum(np.ma.array(data=galaxy_bar, mask=mask))
    
fl_mod[0] = 1
fl_mod_gap[0] = 1   
fl_mod_wei[0] = 1   
    


f, axes = plt.subplots(3,2, figsize=(8,9))

axes[0,0].imshow(galaxy, vmax=1, origin='lower', interpolation='nearest', extent=(-50,50,-50,50))
axes[1,0].imshow(galaxy_hole,vmax=1,  origin='lower', interpolation='nearest', extent=(-50,50,-50,50))
axes[2,0].imshow(galaxy_bar, vmax=1, origin='lower', interpolation='nearest', extent=(-50,50,-50,50))

axes[0,0].set_ylim(-27,27)
axes[0,0].set_xlim(-27,27)

axes[1,0].set_ylim(-27,27)
axes[1,0].set_xlim(-27,27)

axes[2,0].set_ylim(-27,27)
axes[2,0].set_xlim(-27,27)


axes[0,1].plot(rads, fl_mod/max(fl_mod), 'r--')
axes[1,1].plot(rads, fl_mod_gap/max(fl_mod_gap), 'g--')
axes[2,1].plot(rads, fl_mod_wei/max(fl_mod_wei), 'b--')


axes[0,1].set_ylabel('Normalized flux')
axes[1,1].set_ylabel('Normalized flux')
axes[2,1].set_ylabel('Normalized flux')
axes[2,1].set_xlabel('Radius (pixels)')

axes[0,1].set_xlim(5,27)
axes[1,1].set_xlim(5,27)
axes[2,1].set_xlim(5,27)


axes[0,0].set_ylabel(r'$\Delta$Y (pixels)')
axes[1,0].set_ylabel(r'$\Delta$Y (pixels)')
axes[2,0].set_ylabel(r'$\Delta$Y (pixels)')
axes[2,0].set_xlabel(r'$\Delta$X (pixels)')



plt.savefig(ph.MyPATH+'Graphs/Paper_plots/COG_modelling.pdf', bbox_inches='tight')

plt.show()



