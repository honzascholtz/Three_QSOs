
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  2 15:04:55 2021

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
import Tools_path as ph

fsz = gst.graph_format()



Schultze = Table.read(ph.MyPATH+'Catalogues/Schulze19_QSOs.fit_results.fits')

Su_LIR_AGN = Schultze['LIR[agn]']
Su_LIR_SF = Schultze['LIR[gal]']
Su_UV_AGN = Schultze['L2500[agn]']



f, (ax1,ax2) = plt.subplots(2, sharex=False, figsize=(5,10))

ax1.set_xlim(46.3,47.5)
ax2.set_xlim(46.3,47.5)
ax1.set_ylim(46.3,47.5)
ax2.set_ylim(44,47)

ax1.tick_params(direction='in')
ax2.tick_params(direction='in')

ax1.set_xlabel(r'log$_{\rm 10}$ $L_{\rm UV, AGN}$(erg/s)')
ax2.set_xlabel(r'log$_{\rm 10}$ $L_{\rm UV, AGN}$ (erg/s)')
ax1.set_ylabel(r'log$_{\rm 10}$ L$_{\rm IR, AGN}$ (erg/s)')
ax2.set_ylabel(r'log$_{\rm 10}$ L$_{\rm IR, SF}$ (erg/s)')

ax1.plot((46,48), (46,48), linestyle='dotted', color='k')
ax2.plot((46,48), (46,48), linestyle='dotted', color='k')


# =============================================================================
# Extracting and plotting Schukze data points
# =============================================================================
for i in range(len(Su_LIR_AGN)):
    
    El = np.array([Su_LIR_AGN[i,1]-Su_LIR_AGN[i,0] ])
    Eh = np.array([Su_LIR_AGN[i,2]-Su_LIR_AGN[i,1]])
    ers= [El,Eh]
    
    ax1.errorbar(Su_UV_AGN[i,1], Su_LIR_AGN[i,1], yerr=ers, fmt='o', color='k') 
    
    El = np.array([Su_LIR_SF[i,1]-Su_LIR_SF[i,0] ])
    Eh = np.array([Su_LIR_SF[i,2]-Su_LIR_SF[i,1]])
    ers= [El,Eh]
    
    if ers[0]==0:
        ax2.plot(Su_UV_AGN[i,1], ers[1]+0.3, 'ko')
        ax2.plot(Su_UV_AGN[i,1]-0.008, ers[1]+0.21, color='k', marker=arrow, markersize=12)
        
    else:
        ax2.errorbar(Su_UV_AGN[i,1], Su_LIR_SF[i,1], yerr=ers,  fmt='o', color='k')
        

# =============================================================================
#      Extracting and plotting 3 QSOs from this work
        
# =============================================================================
Schultze = Table.read(ph.MyPATH+'Catalogues/QSOs.fit_results.fits')



Su_LIR_AGN = Schultze['LIR[agn]']
Su_LIR_SF = Schultze['LIR[gal]']
Su_UV_AGN = Schultze['L2500[agn]']



for i in range(len(Su_LIR_AGN)):
    
    El = np.array([Su_LIR_AGN[i,1]-Su_LIR_AGN[i,0] ])
    Eh = np.array([Su_LIR_AGN[i,2]-Su_LIR_AGN[i,1]])
    ers= [El,Eh]
    
    ax1.errorbar(Su_UV_AGN[i,1], Su_LIR_AGN[i,1], yerr=ers, fmt='o', color='b') 
    
    El = np.array([Su_LIR_SF[i,1]-Su_LIR_SF[i,0] ])
    Eh = np.array([Su_LIR_SF[i,2]-Su_LIR_SF[i,1]])
    ers= [El,Eh]

    if ers[0]==0:
        ax2.plot(Su_UV_AGN[i,1], ers[1]+0.3, 'bo')
        ax2.plot(Su_UV_AGN[i,1]-0.01, ers[1]+0.21, color='b', marker=arrow, markersize=12)
        
    else:
        ax2.errorbar(Su_UV_AGN[i,1], Su_LIR_SF[i,1], yerr=ers,  fmt='o', color='b')
        

ax1.errorbar((0), (0), fmt='o', color='k', label='Quasars from Schulze et al (2019)')
ax1.errorbar((0), (0), fmt='o', color='b', label='Quasars from this work')

#ax1.set_yticks([46,46.5,47,47.5])

ax1.legend(loc='upper left')

plt.savefig(ph.MyPATH+'Graphs/Paper_plots/SED_test.pdf', bbox_inches='tight')
    

plt.tight_layout()
plt.show()
