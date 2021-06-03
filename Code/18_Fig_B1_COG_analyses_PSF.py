#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 30 15:04:31 2019

@author: jansen
"""

#importing modules
import numpy as np
import matplotlib.pyplot as plt ; plt.ioff()

from astropy.io import fits as pyfits
from astropy import wcs
from astropy.table import Table, join, vstack
from matplotlib.backends.backend_pdf import PdfPages
import pickle
from scipy.optimize import curve_fit

import Graph_setup as gst 
import glob
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
fsz = gst.graph_format(Labelsize=15)

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]

def find_nearest_idx(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

def gauss(x,k,sig):
    mean=0
    expo= -((x-mean)/(np.sqrt(2)*sig))**2    
    return k*e**expo

def straight(x, b):
    y = b*np.ones_like(x)
    return y 

import scipy.optimize as opt

from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.nddata import Cutout2D
import glob

from scipy.misc import imresize
import Tools_plotting as emplot
import Tools_path as ph

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]

def find_nearest_idx(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx



Sample = Table.read(ph.MyPATH+'Four_Quasars.fits')
New = Sample

itere= np.array([2,1,0])


ylims = {'LBQS':( -2.4,2.4)}
ylims['HB89'] = (-2.5,3.2)
ylims['2QZJ'] = (-1,1.25)
ylims['XID_2028']= (-0.06,0.15)


lims_hal={'HB89': np.array([-9,14,-11,12])}
lims_hal['2QZJ'] = np.array([-8,10,-8,10])
lims_hal['LBQS'] = np.array([-10,10,-8,12])
lims_hal['XID_2028'] = np.array([-9,13,-9,13])


column = 0
row = 0

alph = 0.1
axes_width  = 0.20
axes_height = 0.2

im_width  = 0.15
im_height = 0.2

corner_x = 0.03
corner_y = 0.79
    

Fluxes = Table.read(ph.MyPATH+'Catalogues/Target_table_total.fits')

fluxes = np.array([1.3, 1.5, 1.52,0.14])

h = plt.figure(figsize=(15,15))

axes_width = 0.165
corner_x = 0.05
corner_y = 0.88

row= 0 
column = 0


height = 0.1
height_off = 0.12
 
j=0

from scipy.interpolate import interp1d


Halpha_sizes = np.zeros(4)
Halpha_sizes_blr = np.zeros(4)


def smooth(image,sm):
    
    from astropy.convolution import Gaussian2DKernel
    from scipy.signal import convolve as scipy_convolve
    from astropy.convolution import convolve
    
    gauss_kernel = Gaussian2DKernel(sm)

    con_im = convolve(image, gauss_kernel)
    
    #con_im = con_im#*image/image
    
    return con_im  
 


#itere = np.array([3])
for i in itere:
    ID = Sample['ID'][i].strip()
    
    if ID=='HB89':
        IDp='HB8903'
    if ID=='2QZJ':
        IDp='2QZJ00'       
    if ID=='LBQS':
        IDp= 'LBQS01'
    
    if 1==1:
        
        
        
        # Halpha
        
        Image_narrow = pyfits.getdata(ph.MyPATH+'Results_storage/Halpha/'+ID+'_Nearest_spaxel_fit.fits')
        IFU_header = pyfits.getheader(ph.MyPATH+'Results_storage/Halpha/'+ID+'_Nearest_spaxel_fit.fits')   
    
            
        Hal_map = Image_narrow[0]
        
        Hal_map = smooth(Image_narrow[0],1.)
    
        use = np.isnan(Image_narrow[0])
        Hal_map[use] = np.nan
        #Hal_map = smooth(Hal_map,1)
        
        storage = emplot.load_obj('Results_storage/Props/'+ID+'_H')
        storage['Cat_entry'] = New[i]
        storage['X-ray ID'] = ID
        
        shapes = np.shape(Hal_map)
        
            
        IFU_wcs= wcs.WCS(IFU_header).celestial
        
        ax = h.add_axes((corner_x,corner_y-row*height_off, 0.10,height), projection= IFU_wcs)
        ax.tick_params(direction='in')  
        
        cmap = plt.get_cmap('cividis')
        cmap.set_bad(color='black')
        
        if ID=='2QZJ':
            Hal_map = Hal_map*1e18
    
        else:
            Hal_map = Hal_map*1e17
            print ('Normal normalisation')
    
    
        if ID=='2QZJ':
            fl = ax.imshow(Hal_map, cmap=cmap, origin='low',vmax=7., aspect='auto', interpolation='nearest')
    
    
        else:
            fl = ax.imshow(Hal_map, cmap=cmap, origin='low', aspect='auto', interpolation='nearest')
    
        coor = storage['Cube_cent'][1:3].copy()
        coor = np.array(coor, dtype=int)
        lim = lims_hal[ID]
            
        ax.set_xlim(coor[0]+lim[0],coor[0]+lim[1])
        ax.set_ylim(coor[1]+lim[2],coor[1]+lim[3])
            
        ax.set_autoscale_on(False)
        
        dx = lim[1] - lim[0]
        dy = lim[3] - lim[2]
        x = coor[0]+lim[0]
        y = coor[1]+lim[2]
        
        deg_per_pix = IFU_header['CDELT2']
        arc_per_pix = deg_per_pix*3600
        lim_sc = lim*arc_per_pix
        
        ax.text(dx*0.05+x, dy*0.9+y, IDp, color='white')
        
        if (row==2):
            emplot.overide_axes_labels(h,ax,(lim_sc), labelx=1, labely=1,tickin=1, labelsize=15)       
        
        else:
            emplot.overide_axes_labels(h,ax,(lim_sc), labelx=0, labely=1,tickin=1, labelsize=15)
        
        #if row==0:
        #    ax.set_title(r'Narrow H$\alpha$ data', fontsize=15)
            
        if row==0:
            ax.set_title('Data - Multi-fit')
    
    else:
        print( 'No Halpha')
        
# =============================================================================
#  Halpha sub QSO   
# =============================================================================
    
    Image_narrow = np.loadtxt(ph.MyPATH+'Results_storage/Sub_qso_map_hal/'+ID+'.txt')
    IFU_header = pyfits.getheader(ph.MyPATH+'Results_storage/Halpha/'+ID+'_Nearest_spaxel_fit.fits')   
        
    Hal_map = smooth(Image_narrow,1.)
    #Hal_map = Image_narrow
    IFU_wcs= wcs.WCS(IFU_header).celestial
    
    ax = h.add_axes((corner_x+0.13,corner_y-row*height_off, 0.10,height), projection= IFU_wcs)
    
   
    
    if ID=='2QZJ':
        fl = ax.imshow(Hal_map*1e18, cmap=cmap, origin='low',vmin=0.5 , vmax=6., aspect='auto', interpolation='nearest')
        
    else:
        fl = ax.imshow(Hal_map*1e17, cmap=cmap, origin='low',vmin=0.2, aspect='auto', interpolation='nearest')
    
    
    coor = storage['Cube_cent'][1:3].copy()
    coor = np.array(coor, dtype=int)
    lim = lims_hal[ID]
            
    ax.set_xlim(coor[0]+lim[0],coor[0]+lim[1])
    ax.set_ylim(coor[1]+lim[2],coor[1]+lim[3])
            
    ax.set_autoscale_on(False)
        
    dx = lim[1] - lim[0]
    dy = lim[3] - lim[2]
    x = coor[0]+lim[0]
    y = coor[1]+lim[2]
        
    deg_per_pix = IFU_header['CDELT2']
    arc_per_pix = deg_per_pix*3600
    lim_sc = lim*arc_per_pix
        
    ax.text(dx*0.05+x, dy*0.9+y, IDp, color='white')
     
    if (row==2):
        emplot.overide_axes_labels(h,ax,(lim_sc), labelx=1, labely=0,tickin=1, labelsize=15)       
        
    else:
        emplot.overide_axes_labels(h,ax,(lim_sc), labelx=0, labely=0,tickin=1, labelsize=15)
    
    if row==0:
        ax.set_title('Data - BLR-sub')  
# =============================================================================
#     Halpha COG
# =============================================================================
    ax = h.add_axes((corner_x+0.30,corner_y-row*height_off, 0.17,height))
    
    try:
        Res = np.loadtxt(ph.MyPATH+'Results_storage/Growth/'+ID+'Growth_Hal.txt')
    
        ang = Res[0,:]
        
        
        y = Res[1,:]/max(Res[1,:]) #Fluxes['Hal_nar_fl'][i]
        er= Res[3,:]/max(Res[1,:]) #Fluxes['Hal_nar_fl'][i]
        
        if ID=='2QZJ':
            y = Res[1,:]/max(Res[1,:])
            er= Res[3,:]/max(Res[1,:])
            
        
        
        f2 = interp1d(ang, y, kind='cubic')
    
        xpl = np.linspace(min(ang),2,1000)
        ypl = f2(xpl)
    
    
        mid = find_nearest_idx(ypl, 0.5)
        mid = xpl[mid]
        
        Halpha_sizes[i] = mid
        


        yblr = Res[4,:]/Fluxes['Hal_bro_fl'][i]
        
        if ID=='2QZJ':
            yblr = Res[4,:]/Fluxes['Hal_bro_fl'][i]/9
            
        
        
        f2 = interp1d(ang, yblr, kind='cubic')
    
        xpl = np.linspace(min(ang),2,1000)
        ypl = f2(xpl)
    
    
        mid = find_nearest_idx(ypl, 0.5)
        mid = xpl[mid]
        
        Halpha_sizes_blr[i] = mid
        
        
        
        if ID=='XID_2028':
            yblr = Res[4,:]/2.5e-15
            yblr[16:] = 1
        
        nl_hl = ax.plot(Res[0,:],y , label=r'Narrow H$\alpha$', color='red', linewidth=2)    
        ax.fill_between(Res[0,:],y+er, y-er ,color='red', alpha=alph)   
        
        bl_hl = ax.plot(Res[0,:],yblr , label=r'Broad H$\alpha$', color='blue',linestyle='dotted', linewidth=2)    
        
        
    
        #ax.plot(Res[0,:], Res[4,:]/max(Res[4,:]) , label=r'BLR', color='blue', linewidth=2)
    
    
    except:
        print ('No Halpha')
    
    if row==2:
        ax.set_xlabel('radius (arcsec)')
    ax.set_ylabel('Relative flux')  
    
    ax.text(1.0, 0.4, IDp, color='black')
    
    if row==0:
        ax.set_title('Data - COG')  
        
    
# =============================================================================
#     PSF subtraction 
# =============================================================================
    ax = h.add_axes((corner_x+0.52,corner_y-row*height_off, 0.1,height), projection= IFU_wcs)
    
    image = np.loadtxt(ph.MyPATH+'Results_storage/'+ID+'_Halpha_PSF_sub_smt.txt')
    
    if ID=='2QZJ':
        fl = ax.imshow(image*1e18, cmap=cmap, origin='low',vmin=0.5 , vmax=0.83, aspect='auto', interpolation='nearest')
        
    if ID=='LBQS':
        fl = ax.imshow(image*1e17, cmap=cmap, origin='low',vmin=0.2, vmax=0.32, aspect='auto', interpolation='nearest')
    if ID=='HB89':
        fl = ax.imshow(image*1e17, cmap=cmap, origin='low',vmin=0.2, vmax=0.32,  aspect='auto', interpolation='nearest')
               
    ax.set_xlim(coor[0]+lim[0],coor[0]+lim[1])
    ax.set_ylim(coor[1]+lim[2],coor[1]+lim[3])
            
    ax.set_autoscale_on(False)
    
    dx = lim[1] - lim[0]
    dy = lim[3] - lim[2]
    x = coor[0]+lim[0]
    y = coor[1]+lim[2]
        
    deg_per_pix = IFU_header['CDELT2']
    arc_per_pix = deg_per_pix*3600
    lim_sc = lim*arc_per_pix
        
   
    ax.text(dx*0.05+x, dy*0.9+y, IDp, color='white')
    if ID=='HB89':   
        ax.text(dx*0.05+x, dy-3*y, 'Resolved', color='white')
        
    elif ID=='LBQS':   
        ax.text(dx*0.05+x, dy-1*y, 'Unresolved', color='white')
        
    elif ID=='2QZJ':   
        ax.text(dx*0.05+x, y+3, 'Unresolved', color='white')
    
    if (row==2):
        emplot.overide_axes_labels(h,ax,(lim_sc), labelx=1, labely=0,tickin=1, labelsize=15)       
        
    else:
        emplot.overide_axes_labels(h,ax,(lim_sc), labelx=0, labely=0,tickin=1, labelsize=15)
    
    if row==0:
        ax.set_title('PSF residuals') 
    
    
    row+=1
    j+=1



def twoD_Gaussian(ts, amplitude, xo, yo, sigma_x, sigma_y, theta, offset):
    x = ts[0]
    y = ts[1]
    xo = float(xo)
    yo = float(yo)    
    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
    g = offset + amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo) 
                            + c*((y-yo)**2)))
    return g.ravel()


x,y = np.meshgrid(np.arange(100), np.arange(100))


sa = 6
sb = 6
noise = np.random.normal(0, 0.05, 10000).reshape(100,100)
ts=(x,y)
galaxy = twoD_Gaussian(ts, 1., 50,50, sa,sb,0,0).reshape(100,100) #+noise
galaxy_hole = galaxy.copy()
galaxy_bar = galaxy.copy()

x = np.linspace(50-3,50+3, 7)
x = np.array(x, dtype=int)

y = np.linspace(50-3,50+3, 7)
y = np.array(y, dtype=int)

for ix in x:
    for iy in y:
        
        r = np.sqrt((ix-50)**2+(iy-50)**2)
        
        if r<4.:
            galaxy_hole[iy,ix] = galaxy[ix,iy]  - twoD_Gaussian((ix,iy), 1., 50,50, sa,sb,0,0)#.reshape(100,100)


y = np.linspace(50-3,50+3, 7)
y = np.array(y, dtype=int)

x = np.linspace(47,99, 53)
x = np.array(x, dtype=int)

for ix in x:
    for iy in y:        
        galaxy_bar[iy,ix] = galaxy[ix,iy]  - twoD_Gaussian((ix,iy), 1., 50,50, sb,sa,0,0)#.reshape(100,100)

from scipy.misc import imresize


szs = 30
galaxy = imresize(galaxy, (szs, szs), mode='F')
galaxy_hole = imresize(galaxy_hole, (szs, szs), mode='F')
galaxy_bar = imresize(galaxy_bar, (szs, szs), mode='F')

# =============================================================================
# COG
# =============================================================================

N = 20
rads = np.linspace(1,20,N)
fl_mod = np.zeros(N)
fl_mod_gap = np.zeros(N)
fl_mod_wei = np.zeros(N)



mask =  np.zeros((szs,szs)) #np.ma.masked_invalid(galaxy).mask
mask[:,:] = True



for i in range(N):
    mask[:,:] = True
    for ix in range(szs):
        for iy in range(szs):
            dist = np.sqrt((ix- szs/2)**2+ (iy- szs/2)**2)
            if dist< i:
                mask[ix,iy] = False
    
    
    fl_mod[i] = np.ma.sum(np.ma.array(data=galaxy, mask=mask))
    
    fl_mod_gap[i] = np.ma.sum(np.ma.array(data=galaxy_hole, mask=mask))
    
    fl_mod_wei[i] = np.ma.sum(np.ma.array(data=galaxy_bar, mask=mask))
    
fl_mod[0] = 0
fl_mod_gap[0] = 0   
fl_mod_wei[0] = 0   
    
lms = szs/4

# =============================================================================
# Save
# =============================================================================

    
plt.savefig(ph.MyPATH+'Graphs/Paper_plots/Fig4_Halpha_COG.pdf', bbox_inches = 'tight')


Halpha_sizes[0] = -Halpha_sizes_blr[0]
Halpha_sizes[1] = -Halpha_sizes_blr[1]
Halpha_sizes[2] = np.sqrt(Halpha_sizes_blr[2]**2 + Halpha_sizes_blr[2]**2 )


    
plt.show()