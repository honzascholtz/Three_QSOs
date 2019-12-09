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


PATH='/Users/jansen/Google Drive/Astro/'
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

QSO_nz = Table.read(PATH+'Four_Quasars/Catalogues/QSO_Netzer.fits')
QSO_sh = Table.read(PATH+'Four_Quasars/Catalogues/QSO_shen.fits')


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
Kashz_or = Table.read(PATH+'KMOS_SIN/Catalogues/kashz_master_table_v0.3_CIGALE.fits')
Kashz = Kashz_or[np.where(Kashz_or['ZBEST']>1.2)[0]]
Kashz_fortes_or = Table.read(PATH+'KMOS_SIN/Catalogues/kashz.fit_results.fits')
Kashz_fortes = Kashz_fortes_or[np.where(Kashz_fortes_or['Redshift']>1.2)[0]]


Lum = 10**(np.log10(np.array(Kashz['XLUM'], dtype=float))+1)


OIII_k = np.array(Kashz['OIII_LUM_1G'] , dtype=float)
OIII_k[np.where(Kashz['BIC2C_OIII']>2)[0]]= Kashz['OIII_LUM_2GA'][np.where(Kashz['BIC2C_OIII']>2)[0]]+ Kashz['OIII_LUM_2GB'][np.where(Kashz['BIC2C_OIII']>2)[0]]

OIII_k = abs(OIII_k)


# =============================================================================
#  Paper II Sample
# =============================================================================
Sample = Table.read(PATH+'KMOS_SIN/Catalogues/Interesting_Obj_v2.fits')[np.array([0,1,2,4,5,8,10,12])]

Res = Table.read(PATH+'KMOS_SIN/Catalogues/Interesting_Obj_v2_Results_tot.fits')[np.array([0,1,2,4,5,8,10,12])]


Lum_II = 10**(np.log10(np.array(Sample['XLUM'], dtype=float))+1)
Lum_II[-1] = Lum_II[-1]*10

OIII_fl = Res['OIII_nar_fl'] + Res['OIII_bro_fl']

OIII_lum  = np.zeros(8)

for i in range(8):
    
    OIII_lum[i] = luminosity(Sample['ZBEST'][i], OIII_fl[i])



OIII_II = np.array(Sample['OIII_LUM_1G'] , dtype=float)
OIII_II[np.where(Sample['BIC2C_OIII']>2)[0]]= Sample['OIII_LUM_2GA'][np.where(Sample['BIC2C_OIII']>2)[0]]+ Sample['OIII_LUM_2GB'][np.where(Sample['BIC2C_OIII']>2)[0]]

# =============================================================================
#  LOIII vs X-ray lum
# =============================================================================
f,(ax,ax2) = plt.subplots(2, figsize=(8,10))


ax.plot(Lum, OIII_k, '+', color='grey', label='KASHz (z=1.5-2.5)', markersize=8,mew=2, alpha=0.8)
ax.plot(QSO_nz_lum, QSO_nz_lum_oiii, 'g+' , label='QSO Ne+04 (z=1.5-3)', markersize=8,mew=2, alpha=0.8)
ax.plot(QSO_sh_lum, QSO_sh_lum_oiii, 'r+' , label='QSO Sh+16 (z=1.5-3)', markersize=8,mew=2, alpha=0.8)

ax.plot(Lum_II, OIII_II, 'o', color='magenta', label='AGN So+19 (z=1.4-2.6)', markersize=8)


x = np.linspace(40,49,20)
y = x-3.1
ax.plot(10**x, 10**y ,'k--')

x = np.linspace(40,49,20)
y = (x+6.)/1.22

ax.plot(10**x, 10**y ,color='magenta',linestyle='dashed')



xid_fl = 1e-15 
xid_lum = luminosity(1.593,xid_fl)


# Order LBQS, HB89, 2QZJ
Lbol_qso = np.array([QSO_nz_lum[3],QSO_nz_lum[7], QSO_nz_lum[1]])
Loiii_qso = np.array([QSO_nz_lum_oiii[3], QSO_nz_lum_oiii[7], QSO_nz_lum_oiii[1]])

ax.errorbar(1e46,xid_lum, yerr=xid_lum/2,xerr= 5e45 ,color='limegreen', fmt='o', label='XID 2028 (z~1.6)', ms=10)
ax.errorbar(Lbol_qso, Loiii_qso, xerr = Lbol_qso/2 , yerr= Loiii_qso/2, fmt= 'o', color='darkblue', label='QSO - this work (z~2.5)', ms=10)



ax.set_yscale('log')
ax.set_xscale('log')

ax.set_ylim(1e40,9e45)
ax.set_xlim(3e41,2e48)

ax.set_ylabel(r'[OIII] luminosity (ergs s$^{-1}$)')
ax.set_xlabel(r'Bolometric luminosity (ergs s$^{-1}$)')

ax.tick_params(which='both',direction='in')
#ax.legend(loc='upper left',fontsize='large', ncol=2)





# =============================================================================
# FWHM OIIIa
# =============================================================================
W80_k = np.array(Kashz['W80_1G_OIII'] , dtype=float)
W80_k[np.where(Kashz['BIC2C_OIII']>2)[0]]= Kashz['W80_2G_OIII'][np.where(Kashz['BIC2C_OIII']>2)[0]]

W80_II = np.array(Sample['W80_1G_OIII'] , dtype=float)
W80_II[np.where(Sample['BIC2C_OIII']>2)[0]]= Sample['W80_2G_OIII'][np.where(Sample['BIC2C_OIII']>2)[0]]
W80_II[-2] = 700.

FWHMoiii_qso = np.array([QSO_nz_fwhm_oiii[3], QSO_nz_fwhm_oiii[7], QSO_nz_fwhm_oiii[1]])


ax2.plot(OIII_II, W80_II, 'o', color='magenta', markersize=8)

ax2.plot(OIII_k, W80_k, '+', color='grey', markersize=8,mew=2, alpha=0.8)
ax2.plot(QSO_nz_lum_oiii, QSO_nz_fwhm_oiii, 'g+' , markersize=8,mew=2, alpha=0.8)
ax2.plot(QSO_sh_lum_oiii, QSO_sh_fwhm_oiii, 'r+', markersize=8,mew=2, alpha=0.8)

ax2.errorbar(xid_lum, 800., yerr=100,xerr= xid_lum/2 ,color='limegreen', fmt='o', ms=10)
ax2.errorbar(Loiii_qso, FWHMoiii_qso, yerr = FWHMoiii_qso*0.1 , xerr= Loiii_qso/2, fmt= 'o', color='darkblue', ms=10)

ax2.plot(OIII_II, W80_II, 'o', color='magenta', markersize=8)


ax2.tick_params(which='both',direction='in')
ax2.set_ylim(80, 8000)
ax2.set_xlim(1e41,1.5e45)
ax2.set_xscale('log')
ax2.set_yscale('log')

ax2.set_xlabel(r'[OIII] luminosity (ergs s$^{-1}$)')
ax2.set_ylabel(r'[OIII] W80 (km s$^{-1}$)')

ax2.tick_params(which='both',direction='in')
#ax2.legend(loc='upper left',fontsize='large',ncol=2)



legend = f.legend(bbox_to_anchor=(0.19, 0.97), loc=2, borderaxespad=0.,
                   fontsize='12' ,framealpha=1., ncol=2)




plt.savefig(PATH+'Four_Quasars/Graphs/Paper_plots/Fig1_WHM_OIII_Lbol_lum.pdf', bbox='tight')




plt.show()


# =============================================================================
# Spectrum test
# =============================================================================
f, ax = plt.subplots(1)

flux = (1.2, 5.69, 28)
wv = (0.87,3, 2e2)



ax.plot(wv, flux)
ax.set_xscale('log')
ax.set_yscale('log')

