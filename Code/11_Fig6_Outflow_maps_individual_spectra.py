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
import Tools_IFU as IFU
import Tools_plotting as emplot
import Tools_fitting as emfit
import Tools_path as ph

fsz = gst.graph_format()

New= Table.read(ph.MyPATH+'Four_Quasars.fits')



def smooth(image,sm):
    
    from astropy.convolution import Gaussian2DKernel
    from scipy.signal import convolve as scipy_convolve
    from astropy.convolution import convolve
    
    gauss_kernel = Gaussian2DKernel(sm)

    con_im = convolve(image, gauss_kernel)
    
    #con_im = con_im#*image/image
    
    return con_im  



###############################################################################
# Halpha nad ALMA
###############################################################################
k=1

column = 0
row = 0

center = {'LBQS' : np.array([53,47])}
center['HB89'] = np.array([53.7,47.9])
center['2QZJ'] = np.array([53.0, 47.5])

lims={'HB89':  np.array([-10,10,-10,10])}
lims['2QZJ'] = np.array([-10,10,-10,10])
lims['LBQS'] = np.array([-10,10,-10,10])



# =============================================================================
# Halpha, ALMA and and Outflow map
# =============================================================================

f = plt.figure(figsize=(12.1,10))

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
        
    # Halpha    
    
    Maps = pyfits.getdata(ph.MyPATH+'Results_storage/'+ID+'_OIII_mapsind.fits')
    IFU_res = pyfits.getdata(ph.MyPATH+'Results_storage/OIII/'+ID+'_Individual_spaxel_fit.fits')  
    IFU_header = pyfits.getheader(ph.MyPATH+'Results_storage/OIII/'+ID+'_Individual_spaxel_fit.fits')   
    IFU_wcs= wcs.WCS(IFU_header).celestial
    
    # =============================================================================
    # Plotting stuff setup
    # =============================================================================
    if ID=='HB89':    
        
        w80min = 900
        w80max = 1500

        lvls = (1200,1300)
    if ID=='2QZJ':    
        
        w80min = 1600
        w80max = 2000
        
        lvls = (1900,1950)
        
    
    if ID=='LBQS':    
              
        w80min = 1300
        w80max = 2000
        
        lvls = (1800, 1900)
        
        
        
# =============================================================================
#     Plotting images 
# =============================================================================
    ax= f.add_axes((0.035,0.70-(row*0.32), 0.24,0.24), projection= IFU_wcs)
  
    smt = smooth(Maps[2,:,:], 1.5)
    
    use = np.isnan(Maps[2,:,:])
    smt[use] = np.nan
    
    cmap = plt.cm.coolwarm
    cmap.set_bad(color='black')
    
    fl = ax.imshow(smt, vmin=w80min, vmax=w80max, origin='low', cmap=cmap)
    
    ax.contour(smt, levels=lvls, colors='w', linestyles='dashed')
    
    cts=coor = center[ID]
    ax.plot( cts[0], cts[1], 'r*', markersize=12)
    
    lim = lims[ID]
    ax.set_xlim(coor[0]+lim[0],coor[0]+lim[1])
    ax.set_ylim(coor[1]+lim[2],coor[1]+lim[3])
    
    axcbar0 = plt.axes([0.245,0.70-(row*0.32), 0.01,0.24]) 
    axcbar0.tick_params(direction='in')        
    cbar0 = f.colorbar(fl, cax=axcbar0 ,orientation='vertical')
    
# =============================================================================
#     Making it prettier
# =============================================================================
    axcbar0.set_ylabel(r'(km s$^{-1}$)')
    
    deg_per_pix = IFU_header['CDELT2']
    arc_per_pix = deg_per_pix*3600
    
    lim_sc = lim*arc_per_pix
    
    dx = lim[1] - lim[0]
    dy = lim[3] - lim[2]
    x = coor[0]+lim[0]
    y = coor[1]+lim[2]
    
    from matplotlib.patches import Circle
    
    ax.text(dx*0.05+x, dy*0.8+y, IDp, color='white', fontsize=15)
    if ID=='2QZJ':       
        c = Circle((cts[0]-lim[1]*0.6, cts[1]-lim[3]*0.6), 4.5/2, edgecolor='firebrick', facecolor='white', alpha=0.5)
    else:
        c = Circle((cts[0]-lim[1]*0.7, cts[1]-lim[3]*0.7), 4.5/2, edgecolor='firebrick', facecolor='white', alpha=0.5)
     
    ax.add_patch(c)
    
    step=2
    
    ax.plot((coor[0]-step, coor[0]+step), (coor[1]-step, coor[1]-step), color='limegreen', linestyle='solid', linewidth=2)
    ax.plot((coor[0]-step, coor[0]+step), (coor[1]+step, coor[1]+step), color='limegreen', linestyle='solid', linewidth=2)
    
    ax.plot((coor[0]+step, coor[0]+step), (coor[1]-step, coor[1]+step), color='limegreen', linestyle='solid', linewidth=2)
    ax.plot((coor[0]-step, coor[0]-step), (coor[1]-step, coor[1]+step), color='limegreen', linestyle='solid', linewidth=2)
    
    ax.text(coor[0]+step+0.1 , coor[1], '#1', color='white', fontsize=12, fontweight='bold')
    
    if ID=='HB89':
        ax.plot(np.array([coor[0]-step, coor[0]+step])-step*2, np.array([coor[1]-step, coor[1]-step])+step*2, color='orange', linestyle='solid', linewidth=2)
        ax.plot(np.array([coor[0]-step, coor[0]+step])-step*2, np.array([coor[1]+step, coor[1]+step])+step*2, color='orange', linestyle='solid', linewidth=2)
    
        ax.plot(np.array([coor[0]+step, coor[0]+step])-step*2, np.array([coor[1]-step, coor[1]+step])+step*2, color='orange', linestyle='solid', linewidth=2)
        ax.plot(np.array([coor[0]-step, coor[0]-step])-step*2, np.array([coor[1]-step, coor[1]+step])+step*2, color='orange', linestyle='solid', linewidth=2)
    
        ax.text(coor[0]-step+0.1 , coor[1]+2*step, '#2', color='black', fontsize=12, fontweight='bold' )
        
    if ID=='LBQS':
        ax.plot(np.array([coor[0]-step, coor[0]+step])-step*2, np.array([coor[1]-step, coor[1]-step])-step*2, color='orange', linestyle='solid', linewidth=2)
        ax.plot(np.array([coor[0]-step, coor[0]+step])-step*2, np.array([coor[1]+step, coor[1]+step])-step*2, color='orange', linestyle='solid', linewidth=2)
    
        ax.plot(np.array([coor[0]+step, coor[0]+step])-step*2, np.array([coor[1]-step, coor[1]+step])-step*2, color='orange', linestyle='solid', linewidth=2)
        ax.plot(np.array([coor[0]-step, coor[0]-step])-step*2, np.array([coor[1]-step, coor[1]+step])-step*2, color='orange', linestyle='solid', linewidth=2)
        
        ax.text(coor[0]-step+0.1 , coor[1]-2*step, '#2', color='black', fontsize=12, fontweight='bold' )
    
    if ID=='2QZJ':
        ax.plot(np.array([coor[0]-step, coor[0]+step])-step*2, np.array([coor[1]-step, coor[1]-step])+step*2, color='orange', linestyle='solid', linewidth=2)
        ax.plot(np.array([coor[0]-step, coor[0]+step])-step*2, np.array([coor[1]+step, coor[1]+step])+step*2, color='orange', linestyle='solid', linewidth=2)
    
        ax.plot(np.array([coor[0]+step, coor[0]+step])-step*2, np.array([coor[1]-step, coor[1]+step])+step*2, color='orange', linestyle='solid', linewidth=2)
        ax.plot(np.array([coor[0]-step, coor[0]-step])-step*2, np.array([coor[1]-step, coor[1]+step])+step*2, color='orange', linestyle='solid', linewidth=2)
        
        ax.text(coor[0]-step+0.1 , coor[1]+2*step, '#2', color='white', fontsize=12, fontweight='bold' )
        
        
    if row==2:
        emplot.overide_axes_labels(f,ax,(lim_sc), labelx=1, labely=1,tickin=1)

    else:
        emplot.overide_axes_labels(f,ax,(lim_sc), labelx=0, labely=1,tickin=1)
        
    if row==0:
        ax.set_title('W80 map')
    
    row+=1





# =============================================================================
# Their Definition
# =============================================================================
#f = plt.figure(figsize=(8.3,8.3))
print ('# =============================================================================')
print('#Their maps') 
print('# =============================================================================')

itere = np.array([2, 1,0])
row=0
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
    
    Results = pyfits.getdata(ph.MyPATH+'Results_storage/OIII/'+ID+'_Individual_spaxel_fit.fits')
    IFU_header = pyfits.getheader(ph.MyPATH+'Results_storage/OIII/'+ID+'_Individual_spaxel_fit.fits')   
    IFU_wcs= wcs.WCS(IFU_header).celestial
    
# =============================================================================
#     Plotting
# =============================================================================
    ax= f.add_axes((0.35,0.70-(row*0.32), 0.24,0.24), projection= IFU_wcs)
    if ID=='2QZJ':
        z= 2.4063
    if ID=='HB89':
        z = 2.43568286
        
    cnt =  5008.*(1+z)/1e4
    med = (Results[3]-cnt)/cnt*3e5
    
    
    if ID=='LBQS':
        vn = -200
        vx = 300
        
        lvls = (-150,-100,-50)
        
    if ID=='HB89':
        vn = -300  #-700
        vx = -50    #-300
        
        lvls = (-200, -125) #(-550, -500)
        
        if ID=='HB89':
            print ('HB89 Cent ', cnt)
        
        med = med
        med[55,56:63] = -300+300 #-300
        med[54,58:63] = -350+300 #-350
        med[53,60:63] = -400+300 #-400
        
        med[53,60:63] = -400+300 #-400
        med[52,60:63] = -400+300 #-400
        med[51,58:63] = -400+300 #-400
        med[50,60:63] = -400+300 #-400
        #med[50,58:63] = -400
        
    if ID=='2QZJ':
        vn=-550
        vx= -400
        lvls = (-500,-450)
        
    
    smt = smooth(med,1.5)
    
    use = np.isnan(med)
    smt[use] = np.nan
    
   
    fl = ax.imshow(smt,vmin=vn, vmax=vx, origin='low', cmap=cmap)
    
    ax.contour(smt, levels=lvls, colors='white', linewidth=2, linestyle='dashed')
    
    cts=coor = center[ID]
    
    ax.plot( cts[0], cts[1], 'r*', markersize=12)
    
    lim = lims[ID]
    ax.set_xlim(coor[0]+lim[0],coor[0]+lim[1])
    ax.set_ylim(coor[1]+lim[2],coor[1]+lim[3])
  
    
    axcbar0 = plt.axes([0.56,0.70-(row*0.32), 0.01,0.24]) 
    axcbar0.tick_params(direction='in')        
    cbar0 = f.colorbar(fl, cax=axcbar0 ,orientation='vertical')
    
# =============================================================================
#     Making it prettier
# =============================================================================
    axcbar0.set_ylabel(r'(km s$^{-1}$)')
    
    deg_per_pix = IFU_header['CDELT2']
    arc_per_pix = deg_per_pix*3600
    
    lim_sc = lim*arc_per_pix
    
    dx = lim[1] - lim[0]
    dy = lim[3] - lim[2]
    x = coor[0]+lim[0]
    y = coor[1]+lim[2]
    
    from matplotlib.patches import Circle
    
    # PSF
    c = Circle((cts[0]-lim[1]*0.7, cts[1]-lim[3]*0.7), 4.5/2, edgecolor='firebrick', facecolor='white', alpha=0.5)
     
    ax.add_patch(c)
    
    step=2
    
    # Plotting the central rectengular 
    ax.plot((coor[0]-step, coor[0]+step), (coor[1]-step, coor[1]-step), color='limegreen', linestyle='solid', linewidth=2)
    ax.plot((coor[0]-step, coor[0]+step), (coor[1]+step, coor[1]+step), color='limegreen', linestyle='solid', linewidth=2)
    
    ax.plot((coor[0]+step, coor[0]+step), (coor[1]-step, coor[1]+step), color='limegreen', linestyle='solid', linewidth=2)
    ax.plot((coor[0]-step, coor[0]-step), (coor[1]-step, coor[1]+step), color='limegreen', linestyle='solid', linewidth=2)
    
    # Putting label on it 
    ax.text(coor[0]+step+0.1 , coor[1], '#1', color='white', fontsize=12, fontweight='bold' )
    
    # Plotting the off center rectengulars
    if ID=='HB89':
        ax.plot(np.array([coor[0]-step, coor[0]+step])-step*2, np.array([coor[1]-step, coor[1]-step])+step*2, color='orange', linestyle='solid', linewidth=2)
        ax.plot(np.array([coor[0]-step, coor[0]+step])-step*2, np.array([coor[1]+step, coor[1]+step])+step*2, color='orange', linestyle='solid', linewidth=2)
    
        ax.plot(np.array([coor[0]+step, coor[0]+step])-step*2, np.array([coor[1]-step, coor[1]+step])+step*2, color='orange', linestyle='solid', linewidth=2)
        ax.plot(np.array([coor[0]-step, coor[0]-step])-step*2, np.array([coor[1]-step, coor[1]+step])+step*2, color='orange', linestyle='solid', linewidth=2)
    
        ax.text(coor[0]-step+0.1 , coor[1]+2*step, '#2', color='white', fontsize=12, fontweight='bold' )
        
    if ID=='LBQS':
        ax.plot(np.array([coor[0]-step, coor[0]+step])-step*2, np.array([coor[1]-step, coor[1]-step])-step*2, color='orange', linestyle='solid', linewidth=2)
        ax.plot(np.array([coor[0]-step, coor[0]+step])-step*2, np.array([coor[1]+step, coor[1]+step])-step*2, color='orange', linestyle='solid', linewidth=2)
    
        ax.plot(np.array([coor[0]+step, coor[0]+step])-step*2, np.array([coor[1]-step, coor[1]+step])-step*2, color='orange', linestyle='solid', linewidth=2)
        ax.plot(np.array([coor[0]-step, coor[0]-step])-step*2, np.array([coor[1]-step, coor[1]+step])-step*2, color='orange', linestyle='solid', linewidth=2)
        
        ax.text(coor[0]-step+0.1 , coor[1]-2*step, '#2', color='black', fontsize=12, fontweight='bold' )
    
    if ID=='2QZJ':
        ax.plot(np.array([coor[0]-step, coor[0]+step])-step*2, np.array([coor[1]-step, coor[1]-step])+step*2, color='orange', linestyle='solid', linewidth=2)
        ax.plot(np.array([coor[0]-step, coor[0]+step])-step*2, np.array([coor[1]+step, coor[1]+step])+step*2, color='orange', linestyle='solid', linewidth=2)
    
        ax.plot(np.array([coor[0]+step, coor[0]+step])-step*2, np.array([coor[1]-step, coor[1]+step])+step*2, color='orange', linestyle='solid', linewidth=2)
        ax.plot(np.array([coor[0]-step, coor[0]-step])-step*2, np.array([coor[1]-step, coor[1]+step])+step*2, color='orange', linestyle='solid', linewidth=2)
        
        ax.text(coor[0]-step+0.1 , coor[1]+2*step, '#2', color='white', fontsize=12, fontweight='bold' )
    
    
    if row==2:
        emplot.overide_axes_labels(f,ax,(lim_sc), labelx=1, labely=1,tickin=1)

    else:
        emplot.overide_axes_labels(f,ax,(lim_sc), labelx=0, labely=1,tickin=1)
        
    if row==0:
        ax.set_title('Broad component \n offset map')
    
    
    row+=1  
    


row=0
itere = np.array([2,1,0])

print ('# =============================================================================')
print('# Spectra') 
print('# =============================================================================')

#Spectrum lims
spec_up = 5050
spec_do = 4900

txtl = 4910

for i in itere:
    
    print (row)
    ID = New['ID'][i].strip()
      
    z = New['z'][i]    
    if ID=='2QZJ':
        z= 2.4063
    if ID=='HB89':
        z = 2.43568286
    print (i,ID, z)
    
    if ID=='HB89':
        IDp='HB8903'
    if ID=='2QZJ':
        IDp='2QZJ00'       
    if ID=='LBQS':
        IDp= 'LBQS01'
        
    
    # Plotting setup 
    ax_low = f.add_axes((0.69,0.70-(row*0.32), 0.29,0.11))  
    ax_hig = f.add_axes((0.69,0.84-(row*0.32), 0.29,0.11))
    
    ax_low.tick_params(direction='in')
    ax_hig.tick_params(direction='in')
    
    
    row+=1
# =============================================================================
#     Plotting the central spectrum 
# =============================================================================
    # Loading the central squre spectrum 
    Saves = emplot.load_obj('Results_storage/Spectrums/Regions/'+ID+'_OIII_'+str(1)+'_'+str(1))
    
    wave = Saves['wave']/(1+z)*1e4
    Spec = Saves['data']
    
    fit = Saves['total']
    Hbw = Saves['Hbw']
    Hbn = Saves['Hbn']
    
    o3rn = Saves['o3rn']
    o3rw = Saves['o3rw']
    
    o3bn = Saves['o3bn']
    o3bw = Saves['o3bw']
    
    norm=max(fit)
    
    if ID=='2QZJ':
        norm = max(fit)*1.1 # improving scaling for visualization only 
     
    
    # Extracting the whole data including masked
    flux = Spec.data[np.invert(Spec.mask)]
    wv_rst_sc= wave[np.invert(Spec.mask)]
    
    ax_hig.plot(wv_rst_sc, flux/norm, drawstyle='steps-mid', linestyle='solid', color='limegreen')
    
    ax_hig.plot(wave, (o3rn+o3bn)/norm, 'green', linestyle='dashed')
    ax_hig.plot(wave, (o3rw+o3bw)/norm, 'blue', linestyle='dashed')
    
    #ax_hig.plot(wave, (Hbw+Hbn)/max(fit), 'red', linestyle='dotted')
    
    
    ax_hig.set_xlim(spec_do, spec_up)
    ax_hig.set_ylim(-0.1, 1.1)
    
    if ID=='HB89':
         fit = Saves['total']    
         wp = Saves['wave'][np.argmax(fit)]      
         off = (wp-Saves['wave'][np.argmax(o3rw)])/wp*3e5
    else:
        off = Saves['offset']
        wp = Saves['wave'][np.argmax(fit)]  
            
    ax_hig.axvline(wp/(1+z)*1e4, linestyle='--', color='k', alpha=0.3)
    
    ax_hig.text(txtl, 0.85  , r'W80 %.0f km s$^{-1}$' %(Saves['W80'][2]))
    ax_hig.text(txtl, 0.7  , r'Offset %.0f km s$^{-1}$' %(-off))
    
    ax_hig.plot(wave, fit/norm, 'red', linestyle='dotted')
    
    #ax_hig.plot(wave, fit/max(fit), color='red', linestyle='dashed')
    
    ax_low.set_ylabel('                    Normalised flux density')
    
    ax_hig.set_title(IDp)
    
# =============================================================================
#     Plotting the side spectrum 
# =============================================================================
    if ID=='HB89':
        Saves = emplot.load_obj('Results_storage/Spectrums/Regions/'+ID+'_OIII_'+str(0)+'_'+str(0))
    
    if ID=='LBQS':
        Saves = emplot.load_obj('Results_storage/Spectrums/Regions/'+ID+'_OIII_'+str(0)+'_'+str(2))       
        
    if ID=='2QZJ':
        Saves = emplot.load_obj('Results_storage/Spectrums/Regions/'+ID+'_OIII_'+str(0)+'_'+str(0))
    
    wave = Saves['wave']/(1+z)*1e4
    Spec = Saves['data']
    
    fit = Saves['total']
    Hbw = Saves['Hbw']
    Hbn = Saves['Hbn']
    
    o3rn = Saves['o3rn']
    o3rw = Saves['o3rw']
    
    o3bn = Saves['o3bn']
    o3bw = Saves['o3bw']
      
    flux = Spec.data[np.invert(Spec.mask)]
    wv_rst_sc= wave[np.invert(Spec.mask)]
    
    norm = max(fit)
    
    if ID=='2QZJ':
        norm = max(fit)*2
    
    ax_low.axvline(wp/(1+z)*1e4, linestyle='--', color='k', alpha=0.3)
    
    ax_low.plot(wv_rst_sc, flux/norm, drawstyle='steps-mid', linestyle='solid', color='orange')
    
    ax_low.plot(wave, (o3rn+o3bn)/norm, 'green', linestyle='dashed')
    ax_low.plot(wave, (o3rw+o3bw)/norm, 'blue', linestyle='dashed')
    
    ax_low.plot(wave, fit/norm, 'blue', linestyle='dotted')

    
    if ID=='HB89':               
         off = (wp-Saves['wave'][np.argmax(o3rw)])/wp*3e5
    else:
        off = Saves['offset']
    
    ax_low.text(txtl, 0.85  , r'W80 %.0f km s$^{-1}$' %(Saves['W80'][2]))
    ax_low.text(txtl, 0.7  , r'Offset %.0f km s$^{-1}$' %(-off))
    
    
    ax_low.set_xlim(spec_do, spec_up)
    ax_low.set_ylim(-0.1, 1.1)
    
    ax_low.text(5030, 0.8, '#2', color='black')   
    ax_hig.text(5030, 0.8, '#1', color='black')
      
        
ax_low.set_xlabel(r'Restframe wavelength ($\AA$)')
       
    
f.savefig(PATH+'Four_Quasars/Graphs/Paper_plots/Fig5_Outflow_spec_ind.pdf')   


plt.show()