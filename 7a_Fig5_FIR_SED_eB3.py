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
'''
phottable['J_flux']     = np.full(len(maintable),-9999.0)
phottable['J_flux_err'] = np.full(len(maintable),-9999.0)
phottable['H_flux']     = np.full(len(maintable),-9999.0)
phottable['H_flux_err'] = np.full(len(maintable),-9999.0)
phottable['K_flux']     = np.full(len(maintable),-9999.0)
phottable['K_flux_err'] = np.full(len(maintable),-9999.0)
'''
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
# 1.25,1.6,2.2,
band_lams = np.array([3.5,4.6,12.0,22.0,870.0,3000.0, 0.18*1e6]) 


'''
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
'''
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
flux_mor = mor12['mean']
norm_mor = np.interp(6.0,lam_mor,flux_mor)

##  Obtain Mullaney+ 20111 mean and hi luminosity SED 
mullaney11 = Table.read(PATH+'General_data/mullaney_11_seds.tab',
					format='ascii.commented_header',header_start=5)
lam_m11  = mullaney11['Wavelength(um)']
flux_m11_mean  = mullaney11['Mean']/lam_m11
flux_m11_hilum = mullaney11['HiLum']/lam_m11

norm_m11_hilum = np.interp(6.0,lam_m11,flux_m11_hilum)
norm_m11_mean  = np.interp(6.0,lam_m11,flux_m11_mean)

Lyu19 = Table.read(PATH+'General_data/AGN_Lyu+19.txt',
					format='ascii.commented_header',header_start=1,guess=False)

lam_lyu = Lyu19['wavelength(micron)']

lyu_n = Lyu19['NORMAL(jy)']/lam_lyu
lyu_h = Lyu19['HDD(jy)']/lam_lyu
lyu_w = Lyu19['WDD(jy)']/lam_lyu

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

f, axes = plt.subplots(1,1, figsize=(6,5))

i = 0
j = 0


names = ['2QZJ', 'LBQS', 'HB89', 'XID_2028']



mu = np.logspace(3.5, 0.1 , 1000)*1e9
lmb = (c/mu)*1e6 




# Process and plot the photometry normalised to 6 mic 
for iobj in np.array([2]):
    
    print 'Mateos 2012'
    print names[iobj]
    print 'W1-W2', maintable['w1mpro'][iobj]-maintable['w2mpro'][iobj]
    print 'W3-W4', maintable['w3mpro'][iobj]-maintable['w4mpro'][iobj]
    
    ax = axes
    
    z = maintable['z'][iobj]

    fluxes  = np.array(list(phottable[phottable.colnames[0::2]][iobj]))
    efluxes = np.array(list(phottable[phottable.colnames[1::2]][iobj]))

    index, = np.where(fluxes > 0.0)
	
    pfluxes  = (2.9979e14/band_lams[index])*fluxes[index]
    epfluxes = (efluxes[index]/fluxes[index])*pfluxes
    norm = np.interp(6.0,band_lams[index]/(1.0+z),pfluxes)
	
    ax.errorbar(band_lams[index]/(1.0+z),pfluxes/norm,yerr=epfluxes/norm,fmt='ko',ecolor='k')
    ax.errorbar(band_lams[4]/(1.0+z),pfluxes[4]/norm,yerr=epfluxes[4]/norm,fmt='o',color='firebrick',ecolor='firebrick')
    
    # Plot template SEDs using different line colours and styles, normalised to 6 mic
    ax.plot(lam_mor,flux_mor/norm_mor,'r')
    
    ax.plot(lam_lyu,lyu_n/norm_lyu_n,'b')
    ax.plot(lam_lyu,lyu_h/norm_lyu_h,'orange')
    ax.plot(lam_lyu,lyu_w/norm_lyu_w,'yellow')
    
    
    
    xend = lam_mor[-1]
    yend = flux_mor[-1]/norm_mor
    
    ax.plot((xend, 1000), (yend, yend/400) ,'r-')
    
    # SF templates
    ax.plot(sb1_wv, sb1/norm_sb1*(2.9979e14/band_lams[-3])*fluxes[-3]/norm, 'limegreen')
    ax.plot(sb1_wv, sb2/norm_sb2*(2.9979e14/band_lams[-3])*fluxes[-3]/norm, 'darkgreen')
    
    
    print fluxes
    
    #ax.plot(lam_m11,flux_m11_mean/norm_m11_mean,'b--')

    
    
    
    #ax.plot(band_lams[index]/(1.0+z),pfluxes/norm,'k:')
    
    
    # Set the plotting axis ranges
    ax.set_xlim(0.5, 2e5)
    ax.set_ylim(1e-5, 1e1)
    
    ax.set_xscale('log')
    ax.set_yscale('log')
    
    ax.text(1000,1, names[iobj])
    
    index, = np.where(fluxes < 0.0)
    
    pfluxes  = abs((2.9979e14/band_lams[index])*fluxes[index])
    
    
    print index
    ax.plot(band_lams[index]/(1.0+z),pfluxes/norm,marker= u'$\u2193$', markersize=10, color='k')
    
    
    
    flr05 = lmb**(0.5-1)
    flr07 = lmb**(0.7-1)
    
    flr02 = lmb**(-0.2-1)
    flr10 = lmb**(1.-1)
    
    
    
    
    frp = abs((2.9979e14/band_lams[-1])*phottable['VLA_flux'])/norm        
    idx = (np.abs(lmb - (band_lams[-1]/(1.0+z)))).argmin()    
    pfluxes  = abs((2.9979e14/band_lams[-1])*fluxes[-1])
    
    
    
    normc05 = flr05[idx]/(pfluxes/norm)
    normc07 = flr07[idx]/(pfluxes/norm)
    normc02 = flr02[idx]/(pfluxes/norm)
    normc10 = flr10[idx]/(pfluxes/norm)
    
    #ax.fill_between((lmb), flr02/normc02, flr10/normc10 , color='grey', alpha=0.2)    
    ax.fill_between((lmb) ,flr05/normc05, flr07/normc07 , color='magenta', alpha=0.2)
    
    
    flr15 = lmb**(1.5-1)
    
    frp = abs((2.9979e14/band_lams[-2])*phottable['ALMA_Band3_flux'])/norm        
    idx = (np.abs(lmb - (band_lams[-2]/(1.0+z)))).argmin()    
    pfluxes  = abs((2.9979e14/band_lams[-2])*fluxes[-2])
    
    normc15 = flr15[idx]/(pfluxes/norm)
    
    ax.plot((lmb) ,flr15/normc15 , 'b--', alpha=0.7)
    
    ax.tick_params(which='both',direction='in')

    #ax.text(15, 1e-4, '-1')

    ax.set_xlabel(r'Rest Wavelength ($\mu$m)')
        
    
    ax.set_ylabel(r'Normalised luminosity ($\nu L_{\nu}$)')    
        
    j+=1   
    #if j==2:
    #    j=0
    #    i=1




plt.savefig(PATH+'Four_Quasars/Graphs/Paper_plots/HB89_B3_addition.pdf', bbox_inches = 'tight')

Freq = np.array([102.2689,100.4653,  90.25985, 88.52704])
Flux = np.array([5.557, 5.677,  6.708, 6.876])
Error = np.array([ 0.084, 0.064, 0.072, 0.075])



f, ax = plt.subplots(1)


ax.errorbar((Freq), Flux*Freq, yerr=Error*Freq, fmt='go')

ax.set_ylabel(r'$\nu L\nu$ (mJy $\times$ GHz)')
ax.set_xlabel('log(Freq/GHz)')


def fit(freq, k, alpha):    
    y = k*freq**alpha    
    return y
       
    
popt,pcov = curve_fit(fit, Freq, Flux, p0=np.array([100, -1.]))


ax.plot((Freq), (fit(Freq, *popt)*Freq), 'k--')

print popt
plt.savefig(PATH+'Four_Quasars/Graphs/Paper_plots/HB89_B3_slope.pdf', bbox_inches = 'tight')

plt.show()
