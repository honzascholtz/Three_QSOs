import sys
from astropy.table import Table, join, vstack
import numpy as np
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
import astropy.units as u

""" Plotting and examining IR to radio photometry of QSOs from SINFONI/ALMA project
"""


import Graph_setup as gst 
fsz = gst.graph_format()

c= 3e8

plt.close('all')

PATH = '/Users/jansen/Google Drive/Astro/'

# Basis sample is Shimizu+ (2016) catalog 
maintable = Table.read(PATH+'Four_Quasars/Catalogues/4QSO_photometry.fits',format='fits')  # Photometry table compiled by Jan Scholtz
maintable.remove_row(3) # Remove the XID2028 row


phottable = Table()

# Create a table to hold all the photometry
# To expand, add more columns to phottable below
phottable['J_flux']     = np.full(len(maintable),-9999.0)
phottable['J_flux_err'] = np.full(len(maintable),-9999.0)
phottable['H_flux']     = np.full(len(maintable),-9999.0)
phottable['H_flux_err'] = np.full(len(maintable),-9999.0)
phottable['K_flux']     = np.full(len(maintable),-9999.0)
phottable['K_flux_err'] = np.full(len(maintable),-9999.0)
phottable['W1_flux']     = np.full(len(maintable),-9999.0)
phottable['W1_flux_err'] = np.full(len(maintable),-9999.0)
phottable['W2_flux']     = np.full(len(maintable),-9999.0)
phottable['W2_flux_err'] = np.full(len(maintable),-9999.0)
phottable['W3_flux']     = np.full(len(maintable),-9999.0)
phottable['W3_flux_err'] = np.full(len(maintable),-9999.0)
phottable['W4_flux']     = np.full(len(maintable),-9999.0)
phottable['W4_flux_err'] = np.full(len(maintable),-9999.0)
phottable['ALMA_Band7_flux'] = maintable['ALMA Band 7']
phottable['ALMA_Band7_flux_err'] = maintable['ALMA Band 7 er']
phottable['ALMA_Band3_flux'] = maintable['ALMA Band 3']
phottable['ALMA_Band3_flux_err'] = maintable['ALMA Band 3 er']

phottable['VLA_flux'] = maintable['NVSS']
phottable['VLA_flux_err'] = 1.


# observed frame wavelengths of bands in phottable
# Expand as needed, but must have the same order as the bands in phottable
band_lams = np.array([1.25,1.6,2.2,3.5,4.6,12.0,22.0,870.0,3000.0, 0.18*1e6]) 



# Process 2MASS photometry
for iobj in range(len(maintable)):

	if (np.isfinite(maintable['j_m_2mass'][iobj])) & (maintable['j_m_2mass'][iobj] > 0):
		phottable['J_flux'][iobj]      = 1594.0*10.0**(-0.4*maintable['j_m_2mass'][iobj])*1e3 # PS mag + 2MASS tabulated ZP
		phottable['J_flux_err'][iobj]  = phottable['J_flux'][iobj]*maintable['j_msig_2mass'][iobj]

	if (np.isfinite(maintable['h_m_2mass'][iobj])) & (maintable['h_m_2mass'][iobj] > 0):
		phottable['H_flux'][iobj]      = 1024.0*10.0**(-0.4*maintable['h_m_2mass'][iobj])*1e3 # PS mag + 2MASS tabulated ZP
		phottable['H_flux_err'][iobj]  = phottable['H_flux'][iobj]*maintable['h_msig_2mass'][iobj]

	if (np.isfinite(maintable['k_m_2mass'][iobj])) & (maintable['k_m_2mass'][iobj] > 0):
		phottable['K_flux'][iobj]      = 666.7*10.0**(-0.4*maintable['k_m_2mass'][iobj])*1e3 # PS mag + 2MASS tabulated ZP
		phottable['K_flux_err'][iobj]  = phottable['K_flux'][iobj]*maintable['k_msig_2mass'][iobj]

# Process WISE photometry
for iobj in range(len(maintable)):

	if (np.isfinite(maintable['w1mpro'][iobj])) & (maintable['w1sigmpro'][iobj] > 0):
		phottable['W1_flux'][iobj]      = 309.54*10.0**(-0.4*maintable['w1mpro'][iobj])*1e3
		phottable['W1_flux_err'][iobj]  = phottable['W1_flux'][iobj]*maintable['w1sigmpro'][iobj]

	if (np.isfinite(maintable['w2mpro'][iobj])) & (maintable['w2sigmpro'][iobj] > 0):
		phottable['W2_flux'][iobj]      = 171.79*10.0**(-0.4*maintable['w2mpro'][iobj])*1e3
		phottable['W2_flux_err'][iobj]  = phottable['W2_flux'][iobj]*maintable['w2sigmpro'][iobj]

	if (np.isfinite(maintable['w3mpro'][iobj])) & (maintable['w3sigmpro'][iobj] > 0):
		phottable['W3_flux'][iobj]      = 31.67*10.0**(-0.4*maintable['w3mpro'][iobj])*1e3
		phottable['W3_flux_err'][iobj]  = phottable['W3_flux'][iobj]*maintable['w3sigmpro'][iobj]

	if (np.isfinite(maintable['w4mpro'][iobj])) & (maintable['w4sigmpro'][iobj] > 0):
		phottable['W4_flux'][iobj]      = 8.36*10.0**(-0.4*maintable['w4mpro'][iobj])*1e3
		phottable['W4_flux_err'][iobj]  = phottable['W4_flux'][iobj]*maintable['w4sigmpro'][iobj]

##  Obtain Krawczyk QSO template. Not used, since it doesn't extend to the submm
# kraw13 = Table.read('/Users/davidrosario/data/SED_libraries/Krawczyk_2013/apjs468686t2_mrt.txt',format='ascii.cds')
# lam_kraw = 2.9979e14/10**(np.flipud(kraw13['lognu'].data))
# flux_kraw = 10**(np.flipud(kraw13['All']) - 45.0)
# norm_kraw = np.interp(12.0,lam_kraw,flux_kraw)


# =============================================================================
# AGN templates
# =============================================================================
##  Obtain Mor and Netzer 2012 mean SED 
mor12 = Table.read(PATH+'General_data/mor_netzer_mean_and_uncertainty.txt',
					format='ascii.commented_header',header_start=3)
lam_mor = mor12['mic']
flux_mor = mor12['mean']*lam_mor
norm_mor = np.interp(6.0,lam_mor,flux_mor)

##  Obtain Mullaney+ 20111 mean and hi luminosity SED 
mullaney11 = Table.read(PATH+'General_data/mullaney_11_seds.tab',
					format='ascii.commented_header',header_start=5)
lam_m11  = mullaney11['Wavelength(um)']
flux_m11_mean  = mullaney11['Mean']#/lam_m11
flux_m11_hilum = mullaney11['HiLum']#/lam_m11

norm_m11_hilum = np.interp(6.0,lam_m11,flux_m11_hilum)
norm_m11_mean  = np.interp(6.0,lam_m11,flux_m11_mean)

Lyu19 = Table.read(PATH+'General_data/AGN_Lyu+19.txt',
					format='ascii.commented_header',header_start=1,guess=False)

lam_lyu = Lyu19['wavelength(micron)']

lyu_n = Lyu19['NORMAL(jy)']
lyu_h = Lyu19['HDD(jy)']
lyu_w = Lyu19['WDD(jy)']

norm_lyu_n = np.interp(6.0,lam_lyu,lyu_n)
norm_lyu_h = np.interp(6.0,lam_lyu,lyu_h)
norm_lyu_w = np.interp(6.0,lam_lyu,lyu_w)


# =============================================================================
# Star forming tempaltes
# =============================================================================


SED =  Table.read(PATH+'General_data/SED.fits')


sb1_wv = SED['wvl']
sb1 = SED['sb1']/sb1_wv

sb2 = SED['sb2']/sb1_wv

norm_sb1 = np.interp(870./(1+2.5),sb1_wv,sb1)
norm_sb2 = np.interp(870./(1+2.5),sb1_wv,sb2)


phottable['VLA_flux'][:2] = phottable['VLA_flux'][:2]/3

names = ['2QZJ', 'LBQS', 'HB89', 'XID_2028']



mu = np.logspace(3.5, 0.1 , 1000)*1e9
lmb = (c/mu)*1e6 


# =============================================================================
# Templates
# =============================================================================
N = 10
zs = np.linspace(1,4, N)

tps = ['sb1', 'sb2', 'sb3', 'sb4', 'sb5', 'agn']

rats24 = {'dummy':1}


f, ax = plt.subplots(1,figsize=(7,6))

ax.set_ylim(-2.2,3.)
ax.set_xlim(1.4,3)

for j in range(6):
    
    rats = np.zeros(N)
    
    
        
    for i in range(N):
        rats[i] = np.log10(np.interp(870./(1+zs[i]), sb1_wv, SED[tps[j]]) /  np.interp(22./(1+zs[i]), sb1_wv, SED[tps[j]]))
            
    rats24[tps[j]] = rats
    

ax.fill_between(zs, rats24['sb2'], rats24['sb3'], color='green', alpha=0.2,   label='Pure SF - Mul+11')
#ax.fill_between(zs, rats24['agn']-0.2, rats24['agn']+0.2, color='blue', alpha=0.2,   label='Pure AGN - M+11')

rats_hmul = np.zeros(N)
rats_mor = np.zeros(N)

rats_lyun = np.zeros(N)
rats_lyuh = np.zeros(N)
rats_lyuw = np.zeros(N)

for i in range(N):
    rats_hmul[i] = np.log10(np.interp(870./(1+zs[i]), lam_m11,flux_m11_hilum) /  np.interp(22./(1+zs[i]), lam_m11,flux_m11_hilum))    
    rats_mor[i] = np.log10(np.interp(870./(1+zs[i]), lam_mor,flux_mor) /  np.interp(22./(1+zs[i]),lam_mor,flux_mor))
       
    rats_lyun[i] = np.log10(np.interp(870./(1+zs[i]), lam_lyu,lyu_n) /  np.interp(22./(1+zs[i]),lam_lyu,lyu_n))
    rats_lyuh[i] = np.log10(np.interp(870./(1+zs[i]), lam_lyu,lyu_h) /  np.interp(22./(1+zs[i]),lam_lyu,lyu_h))
    rats_lyuw[i] = np.log10(np.interp(870./(1+zs[i]), lam_lyu,lyu_w) /  np.interp(22./(1+zs[i]),lam_lyu,lyu_w))
     
#ax.fill_between(zs, rats_hmul-0.2, rats_hmul+0.2, color='green', alpha=0.2,   label='Pure high AGN - Mul+11')
ers=0.2
ax.fill_between(zs, rats_mor-ers, rats_mor+ers, color='red', alpha=0.2,   label='AGN - Mor+11')

ax.fill_between(zs, rats_lyun-ers, rats_lyun+ers, color='blue', alpha=0.2,   label='AGN - N Lyu+19')
ax.fill_between(zs, rats_lyuh-ers, rats_lyuh+ers, color='orange', alpha=0.2,   label='AGN - H Lyu+19')
ax.fill_between(zs, rats_lyuw-ers, rats_lyuw+ers, color='yellow', alpha=0.2,   label='AGN - W Lyu+19')


#y = -1.19 + 3.623*np.log10(1+zs)
#ax.text(1.1, -0.2, 'strong AGN component')
#ax.plot(zs,y, 'k--')

    
ax.set_ylabel(r'log (F$_{870 \mu \rm m}$/F$_{22 \mu \rm m}$)')
ax.set_xlabel('Redshift')
ax.legend(loc='upper left', ncol=3)

ax.arrow(1.5, -1, 0, 1.5,width=0.007, head_width=0.05, head_length=0.1, fc='k', ec='k')

ax.text(1.55, -0.5, 'Star formation \nproducing 870 \nmicron emission')

f6 = phottable['ALMA_Band7_flux']/ phottable['W4_flux']

ax.errorbar((2.4, 2.3),f6[:2] , yerr= (0.3,0.3), color='green', fmt='o' )
ax.tick_params(direction='in')

f.savefig(PATH+'Four_Quasars/Graphs/Paper_plots/Stanley_plot.pdf',  bbox_inches = 'tight')

plt.show()
