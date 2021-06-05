#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 18 12:04:28 2018

@author: jansen
"""

#importing modules
import numpy as np
import matplotlib.pyplot as plt; plt.ioff()

from astropy.io import fits as pyfits
from astropy import wcs
from astropy.table import Table, join, vstack
from astropy.coordinates import SkyCoord

from matplotlib.backends.backend_pdf import PdfPages
import pickle
from scipy.optimize import curve_fit
import glob


import Graph_setup as gst 
import Tools_IFU as IFU
import Tools_plotting as emplot
import Tools_fitting as emfit
import Tools_path as ph


from scipy.misc import imresize


def smooth(image,sm):
    
    from astropy.convolution import Gaussian2DKernel
    from scipy.signal import convolve as scipy_convolve
    from astropy.convolution import convolve
    
    gauss_kernel = Gaussian2DKernel(sm)

    con_im = convolve(image, gauss_kernel)
    
    #con_im = con_im#*image/image
    
    return con_im  
 
    


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

fsz = gst.graph_format()

New= Table.read(ph.MyPATH+'Four_Quasars.fits')

###############################################################################
# Halpha nad ALMA
###############################################################################
k=1

column = 0
row = 0

# Pixel limits for the images
lims={'HB89': np.array([-9,14,-11,12])}
lims['2QZJ'] = np.array([-20,20,-20,20])
lims['LBQS'] = np.array([-10,10,-8,12])


lims['XID_2028'] = np.array([-9,13,-9,13])


# IFU PSF info. Locations 
IFU_psfs = {'HB89': np.array([-1.1,-0.6,-0.08])/3600}
IFU_psfs['2QZJ'] = np.array([-0.8,-0.7, -0.1 ])/3600
IFU_psfs['LBQS'] = np.array([ 0.0,-0.55, -0.1])/3600

def smooth(image,sm):
    
    from astropy.convolution import Gaussian2DKernel
    from scipy.signal import convolve as scipy_convolve
    from astropy.convolution import convolve
    
    gauss_kernel = Gaussian2DKernel(sm)

    con_im = convolve(image, gauss_kernel)
    
    #con_im = con_im#*image/image
    
    return con_im  

# =============================================================================
# Halpha, ALMA and and Outflow map
# =============================================================================

f = plt.figure(figsize=(6,15.5))

itere = np.array([2, 1,0])

for i in itere:
    ID = New['ID'][i].strip()
      
    z = New['z'][i]    
    
    print (i,ID, z)
    # Halpha    
    if ID=='HB89':
        IDp='HB8903'
    if ID=='2QZJ':
        IDp='2QZJ00'       
    if ID=='LBQS':
        IDp= 'LBQS01'

    Image_narrow = pyfits.getdata(ph.MyPATH+'Results_storage/Halpha/'+ID+'_Nearest_spaxel_fit.fits')
    IFU_header = pyfits.getheader(ph.MyPATH+'Results_storage/Halpha/'+ID+'_Nearest_spaxel_fit.fits')   
    
    # Extracting narrow Hapha flux map and converting it to flux per arcsecond by dividng by (0.1)**2
    Hal_map = Image_narrow[0]
    Hal_map = smooth(Image_narrow[0],1.)/0.01
    
    # The smoothinig spreads the signal to pixel where it wasnt detected. This compensates for it
    use = np.isnan(Image_narrow[0])
    Hal_map[use] = np.nan
    
    storage = emplot.load_obj('Results_storage/Props/'+ID+'_H')
    
    storage['Cat_entry'] = New[i]
    storage['X-ray ID'] = ID
    
    shapes = np.shape(Hal_map)
    
    data = smooth(Hal_map,1)
    data[np.isnan(data)] = 0   
    
    # finding center of the Halpha image
    loc = np.ma.where(data > 0.98*data.max())   
    p_x = np.median(loc[1])
    p_y = np.median(loc[0]) 
        
    IFU_wcs= wcs.WCS(IFU_header).celestial
    
# =============================================================================
#     Plotting
# =============================================================================
    ax= f.add_axes((0.15,0.63-(row*0.29), 0.7,0.27), projection= IFU_wcs)
    
    cmap = plt.get_cmap('cividis')
    cmap.set_bad(color='black')
    
    if ID=='2QZJ':
        fl = ax.imshow(Hal_map*1e16, cmap=cmap, origin='low',vmin=0.5 , vmax=6., aspect='auto', interpolation='nearest')
        
    else:
        fl = ax.imshow(Hal_map*1e15, cmap=cmap, origin='low',vmin=0.2, aspect='auto', interpolation='nearest')
    
    axcbar0 = plt.axes([ 0.85,0.63-(row*0.29) , 0.04, 0.27])  #plt.axes([0.055+ 0.2469,0.397,0.189,0.03])
    axcbar0.tick_params(direction='in')        
    cbar0 = f.colorbar(fl, cax=axcbar0 ,orientation='vertical')
    
    if ID=='2QZJ':  
        axcbar0.set_ylabel(r'SB (10$^{-16}$ erg/s/cm$^{2}$/arcsec$^{2}$)')
    
    else:
        axcbar0.set_ylabel(r'SB (10$^{-15}$ erg/s/cm$^{2}$/arcsec$^{2}$)')
        
        
    coor = storage['Cube_cent'][1:3].copy()
    #coor = np.array(coor, dtype=int)
    lim = lims[ID]
        
    ax.set_xlim(coor[0]+lim[0],coor[0]+lim[1])
    ax.set_ylim(coor[1]+lim[2],coor[1]+lim[3])
        
    ax.set_autoscale_on(False)
    
    ax.scatter(p_x,p_y,marker='o', color='red')
    
    arc = np.round(1./(IFU_header['CDELT2']*3600))
    step = int(arc/4)   
    cent = np.round(storage['Cube_cent'][1:3]).copy()
    cent = np.array(cent, dtype=int)
    cent = cent-0.5
    
    alp=0.4
    '''
    ax.plot( (cent[0]-step, cent[0]-step), (cent[1]+3*step, cent[1]-3*step), 'w-', alpha=alp)
    ax.plot( (cent[0]+step, cent[0]+step), (cent[1]+3*step, cent[1]-3*step), 'w-', alpha=alp)           
    ax.plot( (cent[0]-3*step, cent[0]-3*step), (cent[1]+3*step, cent[1]-3*step), 'w-', alpha=alp)
    ax.plot( (cent[0]+3*step, cent[0]+3*step), (cent[1]+3*step, cent[1]-3*step), 'w-', alpha=alp)
       
    ax.plot( (cent[0]-3*step, cent[0]+3*step), (cent[1]-3*step, cent[1]-3*step), 'w-', alpha=alp)
    ax.plot( (cent[0]-3*step, cent[0]+3*step), (cent[1]+1*step, cent[1]+1*step), 'w-', alpha=alp)     
    ax.plot( (cent[0]-3*step, cent[0]+3*step), (cent[1]-1*step, cent[1]-1*step), 'w-', alpha=alp)
    ax.plot( (cent[0]-3*step, cent[0]+3*step), (cent[1]+3*step, cent[1]+3*step), 'w-', alpha=alp)
    '''
    
    
    if ID=='HB89': # HB89 has strong radio contamination in ALMA band 7 therefore it is as dotted lines
        lnst='dotted'
    else:
        lnst='dashed'
    
    emplot.new_ALMA_contour_plot(storage, ax, both=0, prj = '4QSO', linestyle=lnst)
    
    ax.plot(coor[0], coor[1], 'r*', markersize=10)
    
# =============================================================================
#  Outflow contours  
# =============================================================================
    Maps = pyfits.getdata(ph.MyPATH+'Results_storage/'+ID+'_OIII_mapsind.fits')
    OIII_header = pyfits.getheader(ph.MyPATH+'Results_storage/OIII/'+ID+'_Individual_spaxel_fit.fits')   
    OIII_wcs= wcs.WCS(OIII_header).celestial
    
    smt = smooth(Maps[2,:,:], 1)
    
    use = np.isnan(Maps[2,:,:])
    smt[use] = np.nan
    # =============================================================================
    # Plotting stuff
    # =============================================================================
    if ID=='HB89':    
        lvls = (1200,1300)
        
    if ID=='2QZJ':    
        lvls = (1900,2000)
          
    if ID=='LBQS':    
        lvls = (1800, 1900)
        
        
    ax.contour(smt,transform= ax.get_transform(OIII_wcs), levels=lvls, colors='w', linestyles='dashed')
    
# =============================================================================
#     Outflow Flux contours
# =============================================================================
    Maps = pyfits.getdata(ph.MyPATH+'Results_storage/'+ID+'_OIII_mapsind_flux.fits')
    
    smt = smooth(Maps[1,:,:], 1)
    
    
    use = np.isnan(Maps[1,:,:])
    smt[use] = np.nan
    
    #g,axt = plt.subplots(1)
    #axt.imshow(smt, origin='low')
    
    if ID=='HB89':
        mx= 7e-17
        
        
    if ID=='LBQS':    
        mx= 1.88e-16
          
    if ID=='2QZJ':    
        mx = 8.5e-17
      
    lvls = (0.68*mx, 0.9*mx, 0.95*mx)
        
    ax.contour(smt,transform= ax.get_transform(OIII_wcs), levels=lvls, colors='w', linestyles='solid', linewidths=2)
    
    
# =============================================================================
#     ALMA Position
# =============================================================================
    Img = pyfits.getdata(ph.MyPATH+'ALMA/Final_set/'+ID+'_l_clean.pbcor.fits')[0,0,:,:]
    ALM_header = pyfits.getheader(ph.MyPATH+'ALMA/Final_set/'+ID+'_l_clean.pbcor.fits')  
    
    Ra= New['RA'][i]
    Dec = New['DEC'][i]
    
    ALM_wcs= wcs.WCS(ALM_header).celestial
    
    # Finding the position of the Galaxy in pix scale
    opt_world= np.array([[Ra,Dec]])
    opt_pixcrd = ALM_wcs.wcs_world2pix(opt_world, 0) # WCS transform
    opt_x= (opt_pixcrd[0,0]) # X pixel
    opt_y= (opt_pixcrd[0,1]) # Y pixel
    
    Img[:int(opt_x-30),:] = 0
    Img[int(opt_x+30):,:] = 0
    
    Img[:,:int(opt_x-30)] = 0
    Img[:,int(opt_x+30):] = 0
    
    loc = np.ma.where(Img > 0.98*Img.max())
    
    p_x = np.median(loc[1])
    p_y = np.median(loc[0])
    
    ALM_wcs= wcs.WCS(ALM_header).celestial
    
    alm_world = ALM_wcs.wcs_pix2world(np.array([  [p_x, p_y  ] ]),0)
    alm_world = alm_world[0,:]
    
    ax.scatter(alm_world[0],alm_world[1],transform=ax.get_transform('fk5'), marker='o', color='blue')
       
    
    dx = lim[1] - lim[0]
    dy = lim[3] - lim[2]
    x = coor[0]+lim[0]
    y = coor[1]+lim[2]
    
    ax.text(dx*0.05+x, dy*0.9+y, IDp, color='white', fontsize=15)
    
    deg_per_pix = IFU_header['CDELT2']
    arc_per_pix = deg_per_pix*3600
    
    lim_sc = lim*arc_per_pix
    
    if row==2:
        emplot.overide_axes_labels(f,ax,(lim_sc), labelx=1, labely=1,tickin=1)

    else:
        emplot.overide_axes_labels(f,ax,(lim_sc), labelx=0, labely=1,tickin=1)
    
    row+=1
    
    k+=1


# =============================================================================
#  Legend
# =============================================================================
#ax.plot((-100),(-100), 'ko', label= r'Image Narrow H$\alpha$')
ax.plot((-100),(-100), 'r--', label= 'FIR: star \nformation')
ax.plot((-100),(-100), color='red', linestyle='dotted', label= 'FIR: AGN \nsynchrotron')
ax.plot((-100),(-100), 'w--', label= 'Outflow- velocity')
ax.plot((-100),(-100), 'w-',  label= 'Outflow- SB')

ax.plot((-100),(-100), 'bo', label= r'Centre of H$\alpha$')
ax.plot((-100),(-100), 'ro', label= r'Centre of FIR')
ax.plot((-100),(-100), 'r*', label= 'Quasar \nlocation')


legend = f.legend(bbox_to_anchor=(0.03, 0.995), loc=2, borderaxespad=0.,
                  title=r'Image - Narrow H$\alpha$', fontsize='12' ,framealpha=1., facecolor='black', ncol=3)

plt.setp(legend.get_title(),fontsize='12')

for text in legend.get_texts():
    plt.setp(text, color = 'w')

legend.get_title().set_color('w')

f.savefig(ph.MyPATH+'Graphs/Paper_plots/Fig12_Halpha_ALMA_out_ind.pdf')#), bbox_inches = 'tight')
plt.show()
