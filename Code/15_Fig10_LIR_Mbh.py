#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 20 11:47:25 2020

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

import Tools_path as ph
fsz = gst.graph_format()



Stanley = Table.read(ph.MyPATH+'Catalogues/Stanley_data.fits')
zs = Stanley['z']


Sample = Table.read(ph.MyPATH+'Four_Quasars.fits')
zq = Sample['z'][:-1]

# =============================================================================
# Calculating SMBH masses
# =============================================================================

def MBH_Hbeta(FWHM_bet, L_5100):
   
    Mbh = 6.64+ 2.0*np.log10(FWHM_bet/1000)+ 0.64*(L_5100-44)
    
    return Mbh
    
    
def MBH_Halpha(FWHM_hal, L_Hal):
    
    Mbh = 6.3+ 2.06*np.log10(FWHM_hal/1000)+ 0.55*(L_Hal-42)
    
    
    return Mbh

def luminosity(z,flux):
    from astropy.cosmology import Planck13 as cosmo
    #cosmo = {'omega_M_0':0.3, 'omega_lambda_0':0.7, 'omega_k_0':0.0, 'h':0.72}
    #d_lum= cd.luminosity_distance(z, **cosmo)*(3.08*10**24)
    d_lum = cosmo.luminosity_distance(z)*(3.08*10**24)
    L= 4*pi*(d_lum**2)*flux
    return np.array(L.value, dtype=float) 


def SFR_ms_cal(M,z,const=np.array([0.5,1.5,0.3,0.36,2.5])):
    #Variables
    r= np.log10(1+z)
    m= np.log10(M*10**-9)
    
    #constants    
    m0= const[0]
    a0= const[1]
    a1= const[2]
    m1= const[3]
    a2= const[4]
    
    log10SFR= m - m0 + a0*r - a1*(max(np.array([0., m-m1-a2*r])))**2
    
    return 10**log10SFR

def MBH_Mastar(MBH, z):
    
    Mbh = MBH-8
    
    Mstar = (Mbh+0.52 - np.log10(1+z))/1.12
    
    return Mstar+10

# =============================================================================
# Calculating the Black hole masses from Hbeta and H-alpha
# =============================================================================

L5100 = np.array([3.8e46, 4.9e46, 5e46]) # Carnani 2015
L5100er  = np.array([0.6e46, 0.8e46, 2e46]) # Carnani 2015

FWHM_hal = np.array([6585, 7685, 7488])
FWHM_hal_er = np.array([250, 150, 200])

Hal_fl = np.array([1.8e-15, 8.0e-14, 8.6e-14])
Hal_lum = luminosity(zq, Hal_fl)

Hbet_FWHM = np.array([6295, 6295, 7467])
Hbet_FWHM_er = np.array([200, 200, 200])

Mbh_hal = np.zeros(3)
Mbh_hbe = np.zeros(3)

Mbh_hal_er = np.zeros(3)
Mbh_hbe_er = np.zeros(3)

for i in range(3):
    
    Mbh_Hala = MBH_Halpha( np.random.normal(FWHM_hal[i], FWHM_hal_er[i], 1000), np.random.normal( np.log10(Hal_lum[i]), 0.1, 1000) )
        

    Mbh_Hbea = MBH_Hbeta( np.random.normal(Hbet_FWHM[i], Hbet_FWHM_er[i], 1000), np.random.normal( np.log10(L5100[i]) , 0.3,  1000) )
 
    
    Mbh_hal[i] = np.median(Mbh_Hala)
    Mbh_hbe[i] = np.median(Mbh_Hbea)
    
    
    Mbh_hal_er[i] = np.std(Mbh_Hala)
    Mbh_hbe_er[i] = np.std(Mbh_Hbea)
        
# the masses are consistent with ones estimated by Carniani 2016 - Shemmer 2004
# =============================================================================
# Schulze QSO
# =============================================================================

MBH_Schu = np.array([9.34,9.67, 9.62, 9.48, 9.29, 9.33, 9.49, 9.58, 9.9, 9.59, 9.44, 10.09, 9.39, 9.37, 9.37, 9.41, 9.68, 9.73, 9.54, 9.47 ])
LIR_Schu = np.array([46.06, 45.63, 45.72, 45.1, 44.97, 46.50, 45.5, 45.0, 45.25, 46.01, 45.6, 45.53, 45.17, 45.57, 44.87, 45.67, 45.59, 45, 45.75, 45.73 ])

# =============================================================================
# Making the Plot 
# =============================================================================
Mbh = np.log10(Stanley['M_bh'])
LSF = Stanley['Lir_sf']

LSF_ler =  Stanley['Lir_sf_perr']
LSF_her =   Stanley['Lir_sf']

low = np.where( (zs<0.5) & (zs>0.2) )[0]
mid = np.where( (zs<0.8) & (zs>0.5) )[0]
hig = np.where( (zs<1.5) & (zs>0.8) )[0]
ext = np.where( (zs<3) & (zs>1.5) )[0]

f, ax = plt.subplots(1, figsize=(7.5,6))

# =============================================================================
# Stanley 
Mbh_plt = Mbh[ext]
LSF_her = LSF_her[ext]
LSF_ler = LSF_ler[ext]
LSF = LSF[ext]
srt = Mbh_plt.argsort()

ax.fill_between(Mbh_plt[srt[::-1]], LSF[srt[::-1]]-LSF_ler[srt[::-1]],LSF[srt[::-1]]+LSF_her[srt[::-1]], color='red', alpha=0.3, label = 'Quasars Stanley+17 (1.5<z<2.5) \n- Stacking' )
#ax.errorbar(Mbh[ext], LSF[ext], yerr=[LSF_ler[ext], LSF_her[ext]], fmt='s', color='red', label = 'QSO Stanley+17 (1.5<z<2.5) - Stacking', alpha=0.5, markersize=10)


# =============================================================================
# Schulze
Schulze_SED = Table.read(ph.MyPATH+'Catalogues/Schulze19_QSOs.fit_results.fits')

LIR_schu = np.zeros(len(Schulze_SED))
LIR_schu_do = np.zeros(len(Schulze_SED))
LIR_schu_up = np.zeros(len(Schulze_SED))
LIR_schu_er = np.zeros(len(Schulze_SED))  # this is for the distributions

for i in range(len(Schulze_SED)):
    
    LIR = Schulze_SED['LIR[gal]'][i]
    
    if LIR[1] == 0:
        LIR_schu[i] = LIR[2]
        LIR_schu_er[i] = -99
        
    else:
        LIR_schu[i] = LIR[1]
        LIR_schu_er[i] = 0.3
        
        LIR_schu_do[i] = (10**LIR[1] - 10**LIR[0])
        LIR_schu_up[i] = (10**LIR[2] - 10**LIR[1])
        
        
#LIR_schu_do = 10**LIR_Schu - 10**(LIR_Schu-0.1)
#LIR_schu_up = 10**(LIR_Schu+0.14) - 10**LIR_Schu

ax.errorbar(MBH_Schu, 10**LIR_schu, yerr=[LIR_schu_do, LIR_schu_up], fmt='o',color='green', mfc='white' ,label='Quasars - Schulze+2019', markersize=10)

upper = np.where(LIR_schu_er==-99)[0]
ax.plot(MBH_Schu[upper]-0.008, 10**(LIR_schu[upper]-0.1),marker=arrow, markersize=14,  color='green', linestyle='none')

Schulze_SFR = Table()

Schulze_SFR['SFR'] = LIR_schu-43.41
Schulze_SFR['ESFR'] = 0.3

Schulze_SFR.write(ph.MyPATH+'Catalogues/Schulze_dist.csv', format='csv')



# =============================================================================
#Our QSOs
#QSO_BH = np.array([10.1,10,10.15])

QSO_FIR = 10**np.array([45.16, 45.65 ,44.7+0.48])

QSO_FIR_ler = 10**np.array([44.8, 44.9 , 44.7+0.48])
QSO_FIR_her = 10**np.array([45.9, 46.2 , 44.7+0.48])


QSO_FIR_ler = QSO_FIR - QSO_FIR_ler
QSO_FIR_her = QSO_FIR_her - QSO_FIR

print ('QSO LIR ', QSO_FIR)
print ('QSO LIR ', QSO_FIR_ler)
print ('QSO LIR ', QSO_FIR_her)

print ('QSO SFR ', QSO_FIR*10**-43.41)
print ('QSO SFR ', QSO_FIR_ler*10**-43.41)
print ('QSO SFR ', QSO_FIR_her*10**-43.41)
#Plotting the points
ax.errorbar(Mbh_hbe, QSO_FIR, yerr=[QSO_FIR_ler, QSO_FIR_her], xerr=0.3 , fmt='o', color='darkblue', label='Quasars - This work', markersize=10)
# Plotting the upper limits 
ax.plot(Mbh_hbe[-1]-0.008, QSO_FIR[-1]-QSO_FIR[-1]*0.15,marker = arrow, markersize=14, color='darkblue')

# =============================================================================
# Tracks Main sequence
# =============================================================================
N=20
MBH_plot = np.linspace(8,10.5, N)

Mstar_2 = MBH_Mastar(MBH_plot, 2.0)
Mstar_24 = MBH_Mastar(MBH_plot, 2.4)

LIR_2 = np.zeros(N)
LIR_24 = np.zeros(N)


for i in range(N):

    LIR_2[i] = SFR_ms_cal(10**Mstar_2[i], 2.0)*10**43.41
    LIR_24[i] = SFR_ms_cal(10**Mstar_24[i], 2.4)*10**43.41

ax.plot(MBH_plot, LIR_2, 'g--', label='Main Sequence (z=2.0)', linewidth=2) 
ax.plot(MBH_plot, LIR_24, 'b--',  label='Main Sequence (z=2.4)', linewidth=2) 



# =============================================================================
# Making it pretty
ax.set_yscale('log')
ax.set_xlim(9,10.4)

ax2 = ax.twinx()

yl=ax.get_ylim()

ax2.set_ylim(yl[0]/(10**43.41), yl[1]/(10**43.41))
ax2.set_ylabel(r'SFR (M$_{\odot}$ yr$^{-1}$)')
ax2.set_yscale('log')

ax.legend(loc='upper left', ncol=2)

ax.tick_params(direction='in')
ax.tick_params(direction='in')

ax.set_ylabel(r'L$_{\rm IR, SF}$ (erg/s)', fontsize=20)
ax.set_xlabel(r' log$_{10}$ M$_{\rm BH}$ (M$_{\odot}$)', fontsize=20)

plt.tight_layout()
f.savefig(ph.MyPATH+'Graphs/Paper_plots/Fig8_LIR_MBH.pdf')

plt.show()