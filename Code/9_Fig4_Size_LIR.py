#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 12 14:09:33 2019

@author: jansen
"""

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

PATH='/Users/jansen/Google Drive/Astro/'

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


from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
cosmo = FlatLambdaCDM(H0=72 * u.km / u.s / u.Mpc, Om0=0.3)



# =============================================================================
#  Graph setup
# =============================================================================
f= plt.figure( figsize=(8,8))

ax1 = f.add_axes([0.07,0.08,0.9,0.78])

ax1.set_xscale('log')
ax1.set_ylim(0,9)
ax1.set_xlim(8e44,1e47)
ax1.tick_params(direction='in')

ax1.set_xlabel(r'L$_{\rm FIR,SF}$ (ergs s$^{-1}$)', fontsize=17)
ax1.set_ylabel(r'FIR FWHM (kpc)', fontsize=17)


lfir = np.logspace(44,48,1000)
sfr = lfir*Ken98/Conversion2Chabrier


# =============================================================================
# Schulze data
# =============================================================================
z= 2.
scl = cosmo.kpc_proper_per_arcmin(z).value/60  

LIR_sch = np.array([46, 46.5, 46, 45.6,45.6,45.8,45.7])

Size_sch =    np.array([0.29, 0.26, 0.71, 0.51, 0.33, 0.26, 0.40])*scl
Size_sch_er = np.array([0.07, 0.02, 0.06,0.17, 0.11,0.11,0.11])*scl

ax1.errorbar( 10**LIR_sch, Size_sch, xerr= 0.2*10**LIR_sch, yerr=Size_sch_er, color='darkgreen', linestyle='none' ,fmt='s', markersize=12,  label='Quasars (Schultze+19)', alpha=0.5)



# =============================================================================
# Harrison Data
# =============================================================================
LIR_har = np.array([12.3, 12.6,12.4,12.1,12.3,12.1])
Size_har   = np.array([0.5, 0.2, 0.26, 0.23, 0.17,0.32])
Size_har_er= np.array([0.1, 0.03, 0.04, 0.06, 0.05,0.04])

z_har = np.array([4.1, 1.52, 2.47, 2.39, 1.6, 2.93])




ax1.errorbar(-100,-100,  color='magenta', linestyle='none' ,fmt='s', label='X-ray AGN (Harrison+16)', alpha=0.5)

for i in range(len(LIR_har)):
    
    scl = cosmo.kpc_proper_per_arcmin(z_har[i]).value/60  
    
    print( Size_har[i]/2.35*scl)
    
    xer = [10**LIR_har[i]/2*3.826e33, 10**LIR_har[i]*3.826e33]
    ax1.errorbar(10**LIR_har[i]*3.826e33, Size_har[i]*scl, xerr = 10**LIR_har[i]*3.826e33*0.3, yerr=Size_har_er[i]*scl, markersize=12, color='magenta', linestyle='none' ,fmt='s', alpha=0.5)




# =============================================================================
# Simpson
# ============================================================================
Submm = Table.read(ph.MyPATH+'Catalogues/aS2CLS_sizes_forcmh.fits')



ax1.errorbar(-100, -100,xerr=1, yerr=0 , color='#1f77b4', linestyle='none' ,fmt='s',label='Submm gal (Simpson+15)', alpha=0.5, markersize=12)

for i in range(len(Submm)):
    
    #LIR_er = [(Submm['Lbol_hi'][i]- Submm['Lbol'][i])*3.826e33, (Submm['Lbol'][i]- Submm['Lbol_lo'][i])*3.826e33]
    
    LIR = float(Submm['Lbol'][i])*3.826e33
    
    zs = Submm['Zphot'][i]
    
    
    scl = cosmo.kpc_proper_per_arcmin(zs).value/60  
    
    size = float(Submm['BMAJ'][i])*scl
    size_er = float(Submm['BMAJ_ERR'][i])*scl
    
    ax1.errorbar(LIR,size, xerr=LIR*0.5, yerr=size_er , color='#1f77b4', linestyle='none' ,fmt='o', alpha=0.5, markersize=12)

  
    
# =============================================================================
# Paper II AGN - Scholtz+2020
# =============================================================================
AGN_sizes = np.array([0.9,0.5,2.9,2.7,0.7,1.8])*2.35

AGN_lir = np.array([46.1,46., 45.8,45.7,46.07,46.3])
    
    
ax1.errorbar(10**AGN_lir, AGN_sizes,xerr=0.3*10**AGN_lir, yerr= AGN_sizes*0.2 , color='Firebrick', linestyle='none' ,fmt='s',label='X-ray AGN (Scholtz+20)', alpha=0.3, markersize=12)
   
    
# =============================================================================
# These Quasars
# =============================================================================
QSO_alm = np.array([1.7, 5.5, 3.7,0.6])
QSO_lir=np.array([ 1.4e45 , 4.46e45 ,1.5e45,1.1*9.5e45])


ax1.errorbar(-100, -100,xerr=1, yerr=0 , color='darkblue', linestyle='none' ,fmt='s',markersize=14,label='Quasars - This work')


for i in np.array([0,1]):
    z = Sample['z'][i]
    scl = 1#cosmo.kpc_proper_per_arcmin(z).value/60  
    
    ax1.errorbar(QSO_lir[i], QSO_alm[i]*scl, yerr=QSO_alm[i]*scl*0.1, xerr=0.3*QSO_lir[i], fmt='s',markersize=14,  color='darkblue')
    

for i in np.array([2]):
    z = Sample['z'][i]
    scl = 1# cosmo.kpc_proper_per_arcmin(z).value/60  
    
    ax1.errorbar(QSO_lir[i], QSO_alm[i]*scl, yerr=QSO_alm[i]*scl*0.1, xerr=0, fmt='s',markersize=14,  color='darkblue')
    
    ax1.plot(QSO_lir[i]-0.1*QSO_lir[i],QSO_alm[i]*scl,  marker=u'$\u2190$' , color='darkblue', markersize=13)

for i in np.array([3]):
    z = Sample['z'][i]
    scl = cosmo.kpc_proper_per_arcmin(z).value/60  
    
    ax1.errorbar(QSO_lir[i], QSO_alm[i]*scl, yerr=QSO_alm[i]*scl*0.2, xerr=0.3*QSO_lir[i], fmt='s',markersize=10,  color='limegreen', label='XID 2028 (Brusa+18)')
    
# =============================================================================
# Lmaperti sub
# =============================================================================
    
L_LIR = np.array([229, 184, 44, 85, 362, 48, 384, 108])/Ken98
L_Siz = np.array([0.84, 1, 0.79, 1.4, 1.61, 1.98,2, 2.13 ])*2.35


ax1.errorbar(L_LIR, L_Siz, xerr=(L_LIR*0.5), yerr = L_Siz*0.1,markersize=12,  color='orange',fmt='s', label='X-ray AGN - \n(Lamperti+ in prep)', alpha=0.5)   
    
# =============================================================================
# Fujimoto 2018
# =============================================================================  
Stack_LIR = np.array([1e12, 4e12])*3.8e33
Stack_Siz = np.array([1, 1.05])*2.35

ax1.errorbar(Stack_LIR, Stack_Siz, xerr=(Stack_LIR*0.5), yerr = Stack_Siz*0.1,markersize=12,fmt='o',  color='Firebrick', label='Stacked Galaxies \n(Fujimoto+2018)', alpha=0.5)


# =============================================================================
# Plot end
# =============================================================================
#ax1.legend(loc='upper left', fontsize='large', ncol=2 )

legend = f.legend(bbox_to_anchor=(0.01, 0.99), loc=2, borderaxespad=0.,
                   fontsize='11.5' ,framealpha=1., ncol=3)

plt.savefig(ph.MyPATH+'Graphs/Paper_plots/LIR_size.pdf')#, bbox_inches = 'tight')

plt.show()





