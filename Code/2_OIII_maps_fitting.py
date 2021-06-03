#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 16 14:51:35 2017

@author: jscholtz
"""
from __future__ import print_function
#importing modules
import numpy as np
import matplotlib.pyplot as plt

from astropy.io import fits as pyfits
from astropy import wcs
from astropy.table import Table, join, vstack, Column
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

global_off= 0.2

def luminosity(z,flux,gamma):
    from astropy.cosmology import FlatLambdaCDM
    import astropy.units as u
    cosmo = FlatLambdaCDM(H0=72 * u.km / u.s / u.Mpc, Om0=0.3)
    d_lum = cosmo.luminosity_distance(z)*(3.08*10**24)
    L= 4*pi*(d_lum**2)*flux*(1+z)**(gamma-2)

    return L

fsz = gst.graph_format()


binning = 'Nearest'



import Graph_setup as gst 
import Tools_IFU as IFU
import Tools_plotting as emplot
import Tools_fitting as emfit
import Tools_path as ph


plot_it = 0
#0,4

Sample = Table.read(ph.MyPATH+'Four_Quasars.fits')

Sample['Halpha cube'][0] = '2QZ0028-28'

OIII= True
Spat = False

OIII_size={'HB89': np.array(['0.2', '1'], dtype=str)}
OIII_size['LBQS'] = np.array(['0.2', '1'], dtype=str)
OIII_size['2QZJ'] = np.array(['0.2', '1'], dtype=str)



itere=np.array([2])
for i in itere:
    #ID = New['XID'][i]

    if OIII==True:

        print (ph.MyPATH)
        

        ID = Sample['ID'][i].strip()
        H_file=ID+'_H_fl'
        z= Sample['z'][i]
        OIII_band='Hsin'
        
        if ID=='2QZJ':
            z= 2.4064
            
            

        pdf_plots = PdfPages(ph.MyPATH+ID+'_OIII_tot.pdf')
        
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

        storage_O= IFU.fitting_collapse_OIII(storage_O,z, 1)

        
        Sub_stor=1
        #emplot.Summary_plot(storage_O, 'OIII',z, ID, Sub_stor)
        
        out = storage_O['1D_fit_OIII_mul']
        
        print((out.params['o3rn_center'].value - out.params['o3rw_center'].value)/out.params['o3rn_center'].value*3e5)
        
        pdf_plots.savefig()

        binning = 'Individual'
        storage_O = IFU.Spaxel_fit_mul(storage_O, 'OIII',1, binning,localised=1 )
        
        pdf_plots.close()
       


plt.show()



