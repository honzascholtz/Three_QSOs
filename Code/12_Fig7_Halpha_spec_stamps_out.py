#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  3 19:21:32 2020

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

New= Table.read(PATH+'Four_Quasars/Four_Quasars.fits')

    
import Plotting_tools as emplot

def smooth(image,sm):
    
    from astropy.convolution import Gaussian2DKernel
    from scipy.signal import convolve as scipy_convolve
    from astropy.convolution import convolve
    
    gauss_kernel = Gaussian2DKernel(sm)

    con_im = convolve(image, gauss_kernel)
    
    #con_im = con_im#*image/image
    
    return con_im  


import Graph_setup as gst 
import Tools_IFU as IFU
import Tools_plotting as emplot
import Tools_fitting as emfit
import Tools_path as ph

###############################################################################
# Halpha nad ALMA
###############################################################################
k=1

column = 0
row = 0

lims={'HB89': np.array([-9,14,-11,12])}
lims['2QZJ'] = np.array([-15,15,-15,15])
lims['2QZJ'] = np.array([-20,20,-20,20])

lims['LBQS'] = np.array([-10,10,-8,12])

IFU_psfs = {'HB89': np.array([-1.1,-0.6,-0.08])/3600}
IFU_psfs['2QZJ'] = np.array([-0.8,-0.7, -0.1 ])/3600
IFU_psfs['LBQS'] = np.array([ 0.0,-0.55, -0.1])/3600

BLR = np.array([0.21791792, 0.33963964, 0.32822823, 0.3966967 ])*2.2
BLR[-1] = 0.68   
BLR = BLR/2

PSF=[[0.48, 0.39],
     [0.54, 0.4],
     [0.61, 0.45]]


# =============================================================================
# Halpha, ALMA and and Outflow map
# =============================================================================

f = plt.figure(figsize=(12,10))

# so that we can 


itere = np.array([2,1,0])

print ('# =============================================================================')
print('#Plotting My maps') 
print('# =============================================================================')

for i in itere:
    
    print (row)
    ID = New['ID'][i].strip()
      
    z = New['z'][i]    
    
    print (i,ID, z)
    
    if ID=='HB89':
        IDp='HB8903'
    if ID=='2QZJ':
        IDp='2QZJ00'       
    if ID=='LBQS':
        IDp= 'LBQS01'
        
        
        
    
# =============================================================================
#     Plotting Halpha images
# =============================================================================
    print('# =============================================================================')
    print('# Plotting Halpha maps')
    print('# =============================================================================')
    Image_hal = pyfits.getdata(ph.MyPATH+'Results_storage/Halpha/'+ID+'_Nearest_spaxel_fit.fits')
    Header_hal = pyfits.getheader(ph.MyPATH+'Results_storage/Halpha/'+ID+'_Nearest_spaxel_fit.fits')   
        
    Hal_map = Image_hal[0]
    
    storage = emplot.load_obj('Results_storage/Props/'+ID+'_H')
    
    storage['Cat_entry'] = New[i]
    storage['X-ray ID'] = ID
    
    shapes = np.shape(Hal_map)
    
    data = smooth(Hal_map,1)
    data[np.isnan(data)] = 0   
  
    loc = np.ma.where(data > 0.98*data.max())   

    p_x = np.median(loc[1])
    p_y = np.median(loc[0]) 
        
    Hal_wcs= wcs.WCS(Header_hal).celestial
    
    
    #ax= f.add_axes((0.09+(column*0.48),0.52-(row*0.46), 0.42,0.42), projection= IFU_wcs)
    ax= f.add_axes((0.06,0.70-(row*0.32), 0.2,0.24), projection= Hal_wcs)
    
    cmap = plt.get_cmap('cividis')
    cmap.set_bad(color='black')
    
    
    if ID=='2QZJ':
        ax.imshow(Hal_map, cmap=cmap, origin='low',vmax=7e-18, aspect='auto', interpolation='nearest')
    
    else:
        ax.imshow(Hal_map, cmap=cmap, origin='low', aspect='auto', interpolation='nearest')
    
    coor = storage['Cube_cent'][1:3].copy()
    coor = np.array(coor, dtype=int)
    lim = lims[ID]
        
    ax.set_xlim(coor[0]+lim[0],coor[0]+lim[1])
    ax.set_ylim(coor[1]+lim[2],coor[1]+lim[3])
        
    ax.set_autoscale_on(False)
    
    ax.scatter(p_x,p_y,marker='o', color='red')
    
    arc = np.round(1./(Header_hal['CDELT2']*3600))
    step = int(arc/4)   
    cent = np.round(storage['Cube_cent'][1:3]).copy()
    cent = np.array(cent, dtype=int)
    cent = cent-0.5

# =============================================================================
#  Plotting ALMA
# =============================================================================
    print('# =============================================================================')
    print('# Plotting ALMA Contours')
    print('# =============================================================================')
    
    if ID=='HB89':
        lnst='dotted'
    else:
        lnst='dashed'
        
# =============================================================================
#     Sorting out the plot to make it look pretty
# =============================================================================
    emplot.new_ALMA_contour_plot(storage, ax, both=0, prj = '4QSO', linestyle=lnst)
    
    deg_per_pix = Header_hal['CDELT2']
    arc_per_pix = deg_per_pix*3600
    
    lim_sc = lim*arc_per_pix
    
    emplot.overide_axes_labels(f,ax,(lim_sc), labelx=1, labely=1,tickin=1)

    
# =============================================================================
#  Showing the regional spectra
# =============================================================================
    arc = np.round(1./(Header_hal['CDELT2']*3600))
    step = int(arc/4)   
    
    ax.plot((cent[0]-step, cent[0]+step), (cent[1]-step, cent[1]-step), color='limegreen', linestyle='solid', linewidth=2)
    ax.plot((cent[0]-step, cent[0]+step), (cent[1]+step, cent[1]+step), color='limegreen', linestyle='solid', linewidth=2)
    
    ax.plot((cent[0]+step, cent[0]+step), (cent[1]-step, cent[1]+step), color='limegreen', linestyle='solid', linewidth=2)
    ax.plot((cent[0]-step, cent[0]-step), (cent[1]-step, cent[1]+step), color='limegreen', linestyle='solid', linewidth=2)
    
    if ID=='HB89':
        ax.text(cent[0]-0.5 , cent[1]+step+0.1, '#1', color='white', fontsize=12, fontweight='bold' )   
    else:
        ax.text(cent[0]+step+0.1 , cent[1], '#1', color='white', fontsize=12, fontweight='bold' )
    
    if ID=='HB89':
        ax.plot(np.array([cent[0]-step, cent[0]+step])+step*2, np.array([cent[1]-step, cent[1]-step]), color='orange', linestyle='solid', linewidth=2)
        ax.plot(np.array([cent[0]-step, cent[0]+step])+step*2, np.array([cent[1]+step, cent[1]+step]), color='orange', linestyle='solid', linewidth=2)
    
        ax.plot(np.array([cent[0]+step, cent[0]+step])+step*2, np.array([cent[1]-step, cent[1]+step])+step*0, color='orange', linestyle='solid', linewidth=2)
        ax.plot(np.array([cent[0]-step, cent[0]-step])+step*2, np.array([cent[1]-step, cent[1]+step])+step*0, color='orange', linestyle='solid', linewidth=2)
    
        ax.text(cent[0]+step+1 , cent[1]+step+0.1, '#2', color='white', fontsize=12, fontweight='bold' )
        
    if ID=='LBQS':
        ax.plot(np.array([cent[0]-step, cent[0]+step])-step*2, np.array([cent[1]-step, cent[1]-step])-step*2, color='orange', linestyle='solid', linewidth=2)
        ax.plot(np.array([cent[0]-step, cent[0]+step])-step*2, np.array([cent[1]+step, cent[1]+step])-step*2, color='orange', linestyle='solid', linewidth=2)
    
        ax.plot(np.array([cent[0]+step, cent[0]+step])-step*2, np.array([cent[1]-step, cent[1]+step])-step*2, color='orange', linestyle='solid', linewidth=2)
        ax.plot(np.array([cent[0]-step, cent[0]-step])-step*2, np.array([cent[1]-step, cent[1]+step])-step*2, color='orange', linestyle='solid', linewidth=2)
        
        ax.text(cent[0]-step+0.1 , cent[1]-2*step, '#2', color='white', fontsize=12, fontweight='bold' )
    
    if ID=='2QZJ':
        ax.plot(np.array([cent[0]-step, cent[0]+step])-step*0, np.array([cent[1]-step, cent[1]-step])+step*2, color='orange', linestyle='solid', linewidth=2)
        ax.plot(np.array([cent[0]-step, cent[0]+step])-step*0, np.array([cent[1]+step, cent[1]+step])+step*2, color='orange', linestyle='solid', linewidth=2)
    
        ax.plot(np.array([cent[0]+step, cent[0]+step])-step*0, np.array([cent[1]-step, cent[1]+step])+step*2, color='orange', linestyle='solid', linewidth=2)
        ax.plot(np.array([cent[0]-step, cent[0]-step])-step*0, np.array([cent[1]-step, cent[1]+step])+step*2, color='orange', linestyle='solid', linewidth=2)
        
        ax.text(cent[0]-step+-4 , cent[1]+2*step, '#2', color='white', fontsize=12, fontweight='bold' )
    
    
    dx = lim[1] - lim[0]
    dy = lim[3] - lim[2]
    x = coor[0]+lim[0]
    y = coor[1]+lim[2]
    
    ax.text(dx*0.05+x, dy*0.9+y, IDp, color='white')
    
    
    lim_sc = lim*arc_per_pix
    
    from matplotlib.patches import Circle, Ellipse

    
    c = Circle((coor[0]+lim[1]*0.7, coor[1]-lim[3]*0.3), BLR[i]/arc_per_pix, edgecolor='limegreen', facecolor='limegreen', alpha=0.5)  
    ax.add_patch(c)
    
    PSFs = PSF[0]
    
    c = Ellipse((coor[0]+lim[1]*0.7, coor[1]+lim[3]*0.7), PSFs[0]/arc_per_pix,PSFs[1]/arc_per_pix, 0, edgecolor='firebrick', facecolor='firebrick', alpha=0.5)  
    ax.add_patch(c)  
# =============================================================================
#  setting up the spectra
# =============================================================================
    ax_low = f.add_axes((0.69,0.70-(row*0.32), 0.28,0.24))  
    ax_hig = f.add_axes((0.37,0.70-(row*0.32), 0.28,0.24))
    
    spec_do = 6250
    spec_up = 6900    
    
    
    ax_hig.set_ylabel('Normalised flux density')
    
    ax_low.text(6300, 0.4, '#2', color='black', fontsize=15)   
    ax_hig.text(6300, 0.4, '#1', color='black', fontsize=15)
# =============================================================================
#  Plotting central region         
# =============================================================================

    print('# =============================================================================')
    print('# Plotting Central Halpha Spectra')
    print('# =============================================================================')
    try:
        
        Saves = emplot.load_obj('Results_storage/Spectrums/Regions/'+ID+'_Hal_'+str(1)+'_'+str(1)+'_mod')
        
    except:
        
        Saves = emplot.load_obj('Results_storage/Spectrums/Regions/'+ID+'_Hal_'+str(1)+'_'+str(1)+'_mod')
    
    wv = Saves['wave']
    wave = Saves['wave']/(1+z)*1e4
    Spec = Saves['data']
    
    fit = Saves['total']
    Haw = Saves['Haw']
    Han = Saves['Han']
    
    Nr = Saves['Nr']
    Nb = Saves['Nb']
    
    Hal_fl = Saves['Hal_flx']
    
    flux = Spec.data[np.invert(Spec.mask)]
    wv_rst_sc= wave[np.invert(Spec.mask)]
    
    ax_hig.plot(wv_rst_sc, flux/max(fit), drawstyle='steps-mid', linestyle='solid', color='limegreen')
    
    ax_hig.plot(wave, (Haw)/max(fit), 'magenta', linestyle='dashed')
    ax_hig.plot(wave, (Han)/max(fit), 'orange', linestyle='dashed')
    
    ax_hig.plot(wave, (Nr+Nb)/max(fit), 'green', linestyle='dashed')
    try:
        Hao = Saves['Hao']
        
        ax_hig.plot(wave, (Hao)/max(fit), 'red', linestyle='dashed')
    except:
        lfghs=9999
        
    ax_hig.set_xlim(spec_do, spec_up)
    ax_hig.set_ylim(-0.1, 1.1)
    
    ax_hig.text(6682, 0.7  , 'Narrow H$\\alpha$ flux: \n%.1f $x 10^{-16}$ \nerg/s/cm$^{2}$' %(Hal_fl*1e16))
    
    
    axres = f.add_axes([0.37+0.036, 0.7+0.17-(row*0.32), 0.29/4,0.2/3.5])
               
                
    flux_sc = Spec.data[np.invert(Spec.mask)]
   
    vels = (wv - wv[np.argmax(fit)])/wv[np.argmax(fit)]*3e5
    vels_sc = vels[np.invert(Spec.mask)]
    
    fit_sc = fit[np.invert(Spec.mask)]
    
    
           
    use = np.where((vels_sc<2000) & (vels_sc>-2000))
    
    flux_sc[np.isnan(flux_sc)]=0
    er = np.std(flux_sc-fit_sc)
    
    if ID=='2QZJ':
        er=0.12
    
    axres.plot(vels_sc,flux_sc-fit_sc, drawstyle='steps-mid')
                
    axres.hlines(3*er, -3000,3000, color='k', linestyle='dashed')
    
    axres.hlines(0, -3000,3000, color='r', linestyle='dashed')
                
    axres.text(-750, 3*er,r'3$\sigma$ limit' )
                
                
    axres.set_xlim(-990,990)
    axres.set_ylim(-3*er, 5*er)
    axres.tick_params(direction='in')
    
                
    #axres.set_yticks([])
    axres.set_xlabel('km/s', fontsize=10)

# =============================================================================
#  setting up the spectra
# =============================================================================
    print('# =============================================================================')
    print('# Plotting Extra Halpha Spectra')
    print('# =============================================================================')
    
    
    if ID=='HB89':
        Saves = emplot.load_obj('Results_storage/Spectrums/Regions/'+ID+'_Hal_'+str(2)+'_'+str(1)+'_mod')
    
    if ID=='LBQS':
        Saves = emplot.load_obj('Results_storage/Spectrums/Regions/'+ID+'_Hal_'+str(0)+'_'+str(2)+'_mod')       
        
    if ID=='2QZJ':
        Saves = emplot.load_obj('Results_storage/Spectrums/Regions/'+ID+'_Hal_'+str(1)+'_'+str(0)+'_mod')
    
    wave = Saves['wave']/(1+z)*1e4
    Spec = Saves['data']
    
    fit = Saves['total']
    Haw = Saves['Haw']
    Han = Saves['Han']
    
    Nr = Saves['Nr']
    Nb = Saves['Nb']
    
    Hal_fl = Saves['Hal_flx']
    
    
    
    flux = Spec.data[np.invert(Spec.mask)]
    wv_rst_sc= wave[np.invert(Spec.mask)]
    
    ax_low.plot(wv_rst_sc, flux/max(fit), drawstyle='steps-mid', linestyle='solid', color='orange')
    
    ax_low.plot(wave, (Haw)/max(fit), 'magenta', linestyle='dashed')
    ax_low.plot(wave, (Han)/max(fit), 'orange', linestyle='dashed')
    
    ax_low.plot(wave, (Nr+Nb)/max(fit), 'green', linestyle='dashed')
    
    try:
        Hao = Saves['Hao']
        
        ax_low.plot(wave, (Hao)/max(fit), 'red', linestyle='dashed')
    except:
        lfghs=9999
    
    ax_low.set_xlim(spec_do, spec_up)
    ax_low.set_ylim(-0.1, 1.1)
    
    ax_low.text(6682, 0.7  , 'Narrow H$\\alpha$ flux: \n%.1f $x 10^{-16}$ \nerg/s/cm$^{2}$' %(Hal_fl*1e16))
    
    ax_low.set_xlabel(r'Restframe wavelength ($\AA$)')
    ax_hig.set_xlabel(r'Restframe wavelength ($\AA$)')
    
    axres = f.add_axes([0.69+0.036, 0.7+0.17-(row*0.32), 0.29/4,0.2/3.5])
               
                
    flux_sc = Spec.data[np.invert(Spec.mask)]
   
    vels = (wv - wv[np.argmax(fit)])/wv[np.argmax(fit)]*3e5
    vels_sc = vels[np.invert(Spec.mask)]
    
    fit_sc = fit[np.invert(Spec.mask)]
    
    
           
    use = np.where((vels_sc<2000) & (vels_sc>-2000))
    
    flux_sc[np.isnan(flux_sc)]=0
    er = np.std(flux_sc-fit_sc)
     
    if ID=='2QZJ':
        er=0.05          
    
    axres.plot(vels_sc,flux_sc-fit_sc, drawstyle='steps-mid')
                
    axres.hlines(3*er, -3000,3000, color='k', linestyle='dashed')
    
    axres.hlines(0, -3000,3000, color='r', linestyle='dashed')
                
    axres.text(-750, 3*er,r'3$\sigma$ limit' )
                
                
    axres.set_xlim(-990,990)
    axres.set_ylim(-3*er, 5*er)
                
    axres.tick_params(direction='in')
    
    #axres.tick_params(rotation=45)
             
    axres.set_xlabel('km/s', fontsize=10)
    
    
    row+=1
    


f.savefig(ph.MyPATH+'Graphs/Paper_plots/Fig6_Halpha_maps_spec.pdf')   



plt.show()