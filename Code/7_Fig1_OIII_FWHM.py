#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  1 20:34:11 2019

@author: jansen
"""

#importing modules
import numpy as np
import matplotlib.pyplot as plt;plt.ioff()

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

import Tools_path as ph
fsz = gst.graph_format()

from astropy.cosmology import Planck13 as cosmo
def luminosity(z,flux):
    #cosmo = {'omega_M_0':0.3, 'omega_lambda_0':0.7, 'omega_k_0':0.0, 'h':0.72}
    #d_lum= cd.luminosity_distance(z, **cosmo)*(3.08*10**24)
    d_lum = cosmo.luminosity_distance(z)*(3.08*10**24)
    L= 4*pi*(d_lum**2)*flux
    return float(L.value) 


# =============================================================================
# QSO parent samples
# =============================================================================

QSO_nz = Table.read(ph.MyPATH+'Catalogues/QSO_Netzer.fits')
QSO_sh = Table.read(ph.MyPATH+'Catalogues/QSO_shen.fits')


QSO_nz_lum = 6*(10**np.array(QSO_nz['L5100'], dtype=float))
QSO_nz_lum_oiii = 10**abs(np.array(QSO_nz['col6'], dtype=float))
QSO_nz_fwhm_oiii = np.zeros_like(QSO_nz_lum_oiii)*1.6

for i in range(len(QSO_nz_fwhm_oiii)):    
    try:
        QSO_nz_fwhm_oiii[i] = float(QSO_nz['col8'][i])
    except:
        QSO_nz_fwhm_oiii[i] = np.nan



QSO_sh_lum = 6*(10**np.array(QSO_sh['LOGL5100'], dtype=float))[0]
QSO_sh_lum_oiii = 10**abs(np.array(QSO_sh['LOGL_OIII_5007'], dtype=float))[0]
QSO_sh_fwhm_oiii = np.array(QSO_sh['FWHM_OIII_5007'], dtype=float)[0]*1.6

# =============================================================================
#  KASHz parent sample
# =============================================================================

KASHz_OIII_res = np.loadtxt(ph.MyPATH+'Catalogues/KASHz_OIII.txt')


Lum = KASHz_OIII_res[0,:]
OIII_k =  KASHz_OIII_res[1,:]
W80_k = KASHz_OIII_res[2,:]

# =============================================================================
# SUPER DATA
# =============================================================================
OIII_lum_Super = np.array([43.5, 42.8, 43.16, 42.0, 42.41, 42.24, 43.0, 43.6, 43.3, 42.4, 43.0, 42.9,44.7, 0 , 44.5, 42.8 ,42.7, 42.8 ])


OIII_w80_Super = np.array([2816, 775, 1001, 717, 1044, 1198, 1501, 1015, 1755, 1153, 2142, 1277, 2714, 0, 1457, 678, 1186, 1340])


Lbol_Super = np.array([46.74, 46.8, 46.0,45.44, 46.52,  46.16, 46.82, 46.5, 46.93, 46.0, 46.66, 46.53, 47.91, 47.55, 47.73, 46.5, 46.0, 46.39,])
# =============================================================================
#  Paper II Sample - Scholtz+2020
# =============================================================================
Scholtz_OIII_res = np.loadtxt(ph.MyPATH+'Catalogues/Scholtz_2020_OIII.txt')

Lum_II = Scholtz_OIII_res[0,:]
OIII_II =  Scholtz_OIII_res[1,:]
W80_II = Scholtz_OIII_res[2,:]



# =============================================================================
# SDSS Mullaney 2013
# =============================================================================

SDSS = Table.read(ph.MyPATH+'Catalogues/ALPAKA_W80.fits')

SDSS_XMM = Table.read(ph.MyPATH+'Catalogues/ALPAKA_XMM.fits')
SDSS_XMM_W80 = Table.read(ph.MyPATH+'Catalogues/ALPAKA_W80_XMM.fits')



def luminosity_X(z,flux,gamma):
    #cosmo = {'omega_M_0':0.3, 'omega_lambda_0':0.7, 'omega_k_0':0.0, 'h':0.72}
    #d_lum= cd.luminosity_distance(z, **cosmo)*(3.08*10**24)
    d_lum =( cosmo.luminosity_distance(z)*(3.08*10**24)).value
    d_lum = np.array(d_lum)

    L= 4*pi*(d_lum**2)*flux*(1+z)**(gamma-2)
        
    return L


XFL = SDSS_XMM['EP_3_FLUX'] + SDSS_XMM['EP_4_FLUX'] + SDSS_XMM['EP_5_FLUX']
Xlum = luminosity_X(SDSS_XMM['Z'], XFL , 1.4)

# =============================================================================
#  Harrison 14 data

Lbol_H14 = 10**np.array([45.5,45.0, 45.7, 46.0, 45.6, 46, 45.2, 45.1, 45.3, 45.4, 45.2, 44.3, 45.7, 45.1, 45.3, 45.2])
LOIII_H14 = 10** np.array([42.83, 82.6, 42.80, 43.21, 42.3, 42.82, 42.06, 42.03, 41.93, 42.87, 42.64, 41.95, 42.09, 42.95, 42.72, 42.2])
W80_H14 = np.array([1009, 815, 795, 1449, 1280, 1066, 1285, 778, 1228, 1127, 890, 672, 667, 900, 822, 1184])


# =============================================================================
#  LOIII vs Lbol lum
# =============================================================================
f = plt.figure(figsize=(8,6))

#ax = f.add_axes( (0.1,0.5, 0.85, 0.35 ) )
ax2 = f.add_axes( (0.1,0.1, 0.85, 0.75) )


# Order LBQS, HB89, 2QZJ
xid_fl = 1e-15 
xid_lum = luminosity(1.593,xid_fl)

Lbol_qso = np.array([QSO_nz_lum[3],QSO_nz_lum[7], QSO_nz_lum[1]])
Loiii_qso = np.array([QSO_nz_lum_oiii[3], QSO_nz_lum_oiii[7], QSO_nz_lum_oiii[1]])



# =============================================================================
# SDSS Mullaney 2013 contours
# =============================================================================
binss= list()

binss.append(np.linspace(100, 2000, 15))
binss.append(np.logspace(38,44,20))

hist_sdss, x_sdss,y_sdss  = np.histogram2d(SDSS['OIII_FWHM'], SDSS['OIII_LUM'],bins=binss)
hist_sdss = hist_sdss/np.max(hist_sdss)

# =============================================================================
# FWHM OIIIa
# =============================================================================

FWHMoiii_qso = np.array([1260,1980,1900]) # number from our OIII fitting (see OIII_regions script)

#ax2.errorbar(xid_lum, 800., yerr=100,xerr= xid_lum/2 ,color='limegreen', fmt='o', ms=10)
l1 = ax2.errorbar(Loiii_qso, FWHMoiii_qso, yerr = FWHMoiii_qso*0.1 , xerr= Loiii_qso/2, fmt= 'o', color='darkblue', ms=10, label='Quasars - this work (z~2.5)')

l3 = ax2.plot(QSO_nz_lum_oiii, QSO_nz_fwhm_oiii, 'g+' , markersize=8,mew=2, alpha=0.8, label='Quasars - Parent sample (z=1.5-3)')
ax2.plot(QSO_sh_lum_oiii, QSO_sh_fwhm_oiii, 'g+', markersize=8,mew=2, alpha=0.8)


ax2.contour(y_sdss[:-1],x_sdss[:-1], 1-hist_sdss, levels= (0.1, 0.68,0.8, 0.9, 0.95), colors= 'orange', linewidths=2, alpha=0.7, label='AGN- SDSS (z<0.4)')
l4 = ax2.plot(1,1, color= 'orange', linewidth=2, alpha=0.7, label='AGN - SDSS (z<0.4)')

l7 = ax2.plot(OIII_II, W80_II, 'o', color='magenta', markersize=8, label='X-AGN Sch+20 (z=1.4-2.6)')

l6 = ax2.plot(OIII_k, W80_k, '+', color='grey', markersize=8,mew=2, alpha=0.8, label='X-AGN KASHz (z=1.5-2.5)')

l8 = ax2.plot(10**OIII_lum_Super,OIII_w80_Super, 'ro' , label='X-AGN SUPER (z$\sim 2$)')


# =============================================================================
#Harrison point
#ax2.plot(LOIII_H14, W80_H14, 'o', color='firebrick', linestyle='none')


ax2.tick_params(which='both',direction='in')
ax2.set_ylim(100, 4000)
ax2.set_xlim(2e38,1.5e45)
ax2.set_xscale('log')
ax2.set_yscale('log')

ax2.set_xlabel(r'[OIII] luminosity (ergs s$^{-1}$)')
ax2.set_ylabel(r'[OIII] W80 (km s$^{-1}$)')

ax2.tick_params(which='both',direction='in')
#ax2.legend(loc='upper left',fontsize='large',ncol=2)


handles, labels = ax2.get_legend_handles_labels()
order = [5, 0,1,2,4,3]

print(labels)

legend = f.legend([handles[idx] for idx in order],[labels[idx] for idx in order], bbox_to_anchor=(0.12, 0.99), loc=2, borderaxespad=0.,
                   fontsize='12' ,framealpha=2., ncol=2)

#ax2.legend()
plt.savefig(ph.MyPATH+'Graphs/Paper_plots/Fig1_FWHM_OIII.pdf', bbox='tight')



z=2.5
d = cosmo.luminosity_distance(z)*(3.08*10**22)

f14 = 29.4 /1000 * 1e-26
L = 4*np.pi*d*d * f14 *(1+z)**-(1-0.7)
print (L)

plt.show()

