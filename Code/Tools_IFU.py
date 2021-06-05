#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  3 16:36:26 2017

@author: jansen
"""

#importing modules
import numpy as np
import matplotlib.pyplot as plt

from astropy.io import fits
from astropy.wcs import wcs
from astropy.nddata import Cutout2D

import Tools_fitting as emfit
import Tools_path as ph
#import Tools_plotting as emplot
from matplotlib.backends.backend_pdf import PdfPages




nan= float('nan')

pi= np.pi
e= np.e

c= 3.*10**8
h= 6.62*10**-34
k= 1.38*10**-23

Ken98= (4.5*10**-44)
Conversion2Chabrier=1.7 # Also Madau
Calzetti12= 2.8*10**-44
arrow = u'$\u2193$' 




def conf(aray):
    sorted_array= np.array(sorted(aray))
    leng= (float(len(aray))/100)*16
    leng= int(leng)
    
    
    hgh = sorted_array[-leng]
    low = sorted_array[leng]
    
    return low, hgh


def load_cube(Full_path, z, ID, flag):
    '''Loads the cube and extracts following information:
        header,
        dim - dimensions
        flux
        noise
    '''
    
    filemarker = fits.open(Full_path)
    
    print (Full_path)
    if flag=='KMOS':
        
        header = filemarker[1].header # FITS header in HDU 1
        flux_temp  = filemarker[1].data/1.0e-13   
                         
        filemarker.close()  # FITS HDU file marker closed
    
    elif flag=='Sinfoni':
        header = filemarker[0].header # FITS header in HDU 1
        flux_temp  = filemarker[0].data/1.0e-13*1e4
                         
        filemarker.close()  # FITS HDU file marker closed
    
    else:
        print ('Instrument Flag is not understood!')
    
    flux  = np.ma.masked_invalid(flux_temp)   #  deal with NaN
   
    # Number of spatial pixels
    n_xpixels = header['NAXIS1']
    n_ypixels = header['NAXIS2']
    dim = [n_ypixels, n_xpixels]

    #  Number of spectral pixels
    n_spixels = header['NAXIS3']
    dim = [n_ypixels, n_xpixels, n_spixels]
    
    
     # Fixing issue with keywords in the header   
    try:
        x = header['CDELT1']
    except:
        header['CDELT1'] = header['CD1_1']
    
    
    try:
        x = header['CDELT2']
    except:
        header['CDELT2'] = header['CD2_2']
        
    
    try:
        x = header['CDELT3']
    except:
        header['CDELT3'] = header['CD3_3']
    
    wave = header['CRVAL3'] + (np.arange(n_spixels) - (header['CRPIX3'] - 1.0))*header['CDELT3']
    
    storage = {'dim': dim}
    storage['obs_wave'] = wave
    storage['flux'] = flux
    storage['z_guess'] = z
    storage['X-ray ID'] = ID
    storage['Telescope_flag'] = flag
    
    try:
        deg_per_pix_x = abs(header['CDELT1'])
        
    except:
        deg_per_pix_x = header['CDELT1']
        
        
    arc_per_pix_x = 1.*deg_per_pix_x*3600
    Xpix = header['NAXIS1']
    Xph = Xpix*arc_per_pix_x
    
    try:
        deg_per_pix_y = abs(header['CDELT2'])
    
    except:
        deg_per_pix_y = header['CDELT2']
        header['CD2_2'] = header['CDELT2']
        
    arc_per_pix_y = deg_per_pix_y*3600
    Ypix = header['NAXIS2']
    Yph = Ypix*arc_per_pix_y
    
    storage['header'] =  header
    storage['Phys_size'] = np.array([Xph, Yph])
    
    return storage


def add_res(storage, line_cat):
    
    storage['Cat_entry'] = line_cat
    
    return storage


def mask_emission(storage, z):
    '''This function masks out all the OIII and HBeta emission
    '''
    OIIIa=  501./1e3*(1+z)
    OIIIb=  496./1e3*(1+z)
    Hbeta=  485./1e3*(1+z)
    width = 0.006

    mask =  storage['flux'].mask.copy()
    
    wave= storage['obs_wave']


    
    OIIIa_loc = np.where((wave<OIIIa+width)&(wave>OIIIa-width))[0]
    OIIIb_loc = np.where((wave<OIIIb+width)&(wave>OIIIb-width))[0]
    Hbeta_loc = np.where((wave<Hbeta+width)&(wave>Hbeta-width))[0]

    mask[OIIIa_loc,:,:] = True
    mask[OIIIb_loc,:,:] = True
    mask[Hbeta_loc,:,:] = True
    
    storage['em_line_mask'] = mask
    
    return storage



def mask_sky(storage,sig, mode=0):
        
    bins = np.linspace(0,2048,5)
    bins= np.array(bins, dtype=int)
    
    flux = np.ma.array(data=storage['flux'].data.copy(), mask= storage['em_line_mask'].copy())    
    mask =  storage['em_line_mask'].copy()  
    dim = storage['dim']
    wave = storage['obs_wave']
    x=0
    
    if mode=='Hsin':
        use = np.where((wave>1.81)&(wave<1.46))[0]
        mask[use,:,:] = True 
    
    for i in range(dim[0]):
        for j in range(dim[1]):
            stds = np.ma.std(flux[:,i,j])
            
            try:
                p_coeff = np.polyfit( wave, flux[:,i,j], 5 ) 
                f2 = np.poly1d( p_coeff ) 
                y = f2(wave)
                     
                sky = np.where((flux[:,i,j]< (y-stds*sig)) | (flux[:,i,j]> (y+stds*sig)))[0]
                #print sky
                mask[sky,i,j] = True
                
            except:
                x+=1
                #print x
                #print i,j
        
            
    storage['sky_line_mask-em'] = mask
    
    flux = np.ma.array(data=flux.data, mask=mask)
    
    noise_spax = np.ma.std(flux, axis=(0))
            
    storage['spax_noise'] = noise_spax
    
   
    return storage


def unmask_em(storage, z):
    ''' This unmasks the the Emission lines
    '''

    OIIIa=  501./1e3*(1+z)
    OIIIb=  496./1e3*(1+z)
    Hbeta=  485./1e3*(1+z)
    width = 0.006

    mask = storage['sky_line_mask-em']
    
    wave= storage['obs_wave']

    OIIIa_loc = np.where((wave<OIIIa+width)&(wave>OIIIa-width))[0]
    OIIIb_loc = np.where((wave<OIIIb+width)&(wave>OIIIb-width))[0]
    Hbeta_loc = np.where((wave<Hbeta+width)&(wave>Hbeta-width))[0]

    mask[OIIIa_loc,:,:] = False
    mask[OIIIb_loc,:,:] = False
    mask[Hbeta_loc,:,:] = False
    
    storage['sky_clipped+em'] = mask
    
    return storage 


def collapse_white(storage, plot):
    flux = np.ma.array(data=storage['flux'].data, mask=storage['sky_line_mask-em']) 
    
    median = np.ma.median(flux, axis=(0))  
    
    ID = storage['X-ray ID']
    if ID== 'ALESS_75':
        wave = storage['obs_wave'].copy()
        
        use = np.where((wave < 1.7853) &(wave > 1.75056))[0]
        
        use = np.append(use, np.where((wave < 2.34) &(wave > 2.31715))[0] )
        
        median = np.ma.median(flux[use,:,:], axis=(0))  
        
   
    storage['Median_stack_white'] = median
    
  
    if plot==1:
        plt.figure()
        plt.imshow(median,  origin='low')
        plt.colorbar()
        
    return storage


def twoD_Gaussian(dm, amplitude, xo, yo, sigma_x, sigma_y, theta, offset): 
    x, y = dm
    xo = float(xo)
    yo = float(yo)    
    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
    g = offset + amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo) 
                            + c*((y-yo)**2)))
    return g.ravel()



def find_center(storage, im, plot, extra_mask=0, manual=np.array([0])):
    '''
    Input: 
        Storage, 
        image name to be loaded
        Plot it?
        extra mask of the weird bright features
        manual - If manual is on then It will just return the manual inputs. Use if there is no continuum detected
    '''
    
    shapes = storage['dim']
    data = storage[im].copy()
    
    # Create the x inputs for the curve_fit
    x = np.linspace(0, shapes[1]-1, shapes[1])
    y = np.linspace(0, shapes[0]-1, shapes[0])
    x, y = np.meshgrid(x, y)   
    
    import scipy.optimize as opt
    
    if len(manual)==1:
        
        # If there is no extra mask -create a dummy extra mask
        if len(np.shape(extra_mask))<1.:
            extra_mask = storage['flux'].mask.copy()
            extra_mask = extra_mask[0,:,:]
            extra_mask[:,:] = False
    
        #############
        # Masking the edges based on the pixel scale. 
        edges = storage['flux'].mask.copy() # Mask to collapse
        edges = edges[0,:,:]
        edges[:,:] = False
        try:
            pixel_scale = 1./(storage['header']['CD2_2']*3600)
        
        except:
            pixel_scale = 1./(storage['header']['CDELT2']*3600)
            
    
        if pixel_scale < 7:    
            edges[:,0] = True
            edges[:,1] = True
            edges[:, -1] = True
            edges[:, -2] = True
            edges[1,:] = True
            edges[0,:] = True
            edges[-1,:] = True
            edges[-2,:] = True
            
            print ('Masking edges based on 0.2 scale')
            
        else:
            edges[:,:5] = True
            edges[:, -6:] = True        
            edges[:5,:] = True        
            edges[-6:,:] = True
            
            print ('Masking edges based on 0.1 scale')
        
        # Combining the edge and extra mask 
        comb = np.logical_or(extra_mask, edges)
        
        # Setting all other nan values to 0 to avoid any troubles
        data.data[np.isnan(data.data)] = 0
        data= data.data
        
        # combining the data and masks
        masked_im = np.ma.array(data= data, mask=comb)
        
        # Finding the center of the contiunuum 
        loc = np.ma.where(masked_im == masked_im.max())
        print ('location of the peak on the continuum',loc)
        
        #plt.figure()
        #plt.title('Find center - plot')
        #plt.imshow(masked_im, origin='low')
        
        # Setting the            
        initial_guess = (data[loc[1][0], loc[0][0]],loc[1][0],loc[0][0],1,1,0,0)
        
        print ('Initial guesses', initial_guess)
        dm = (x,y)
        popt, pcov = opt.curve_fit(twoD_Gaussian, dm, data.ravel(), p0=initial_guess)
        
        er = np.sqrt(np.diag(pcov))
    
        print ('Cont loc ', popt[1:3])
        print ('Cont loc er', er[1:3])
        storage[im+'_Center_data'] = popt 
        
        plt.show()
        
    else:
        manual = np.append(data[int(manual[0]), int(manual[1])], manual)
        manual = np.append(manual, np.array([2.,2.,0.5,0. ]))
        storage[im+'_Center_data'] = manual
    
    return storage

def choose_pixels(storage, plot, rad= 0.6, flg=1):
    ''' Choosing the pixels that will collapse into the 1D spectrum. Also this mask is used later
    '''
    center =  storage['Median_stack_white_Center_data'][1:3].copy()
    shapes = storage['dim']
    
    print ('Center of cont', center)
    
    print ('Extracting spectrum from diameter', rad*2, 'arcseconds')
    
    # Creating a mask for all spaxels. 
    mask_catch = storage['flux'].mask.copy()
    mask_catch[:,:,:] = True
    header  = storage['header']
    #arc = np.round(1./(header['CD2_2']*3600))
    arc = np.round(1./(header['CDELT2']*3600))
    print ('radius ', arc*rad)
    
    # This choose spaxel within certain radius. Then sets it to False since we dont mask those pixels
    for ix in range(shapes[0]):
        for iy in range(shapes[1]):
            dist = np.sqrt((ix- center[1])**2+ (iy- center[0])**2)
            if dist< arc*rad:
                mask_catch[:,ix,iy] = False
    
    # 587 have special extended stuff that I dont want to mask
    if flg=='K587':
        mask_catch[:,11:29 ,3:36 ] = False
        print ('extracting flux K band XID 587')
    
    elif flg=='H587':
        mask_catch[:,12:29 ,2:36 ] = False
        print ('extracting flux H band XID 587')
    
    elif flg=='K751':
        mask_catch[:,3:31 ,5:34 ] = False
        print ('extracting flux K band XID 751')
    
    elif flg=='K587sin':
        mask_catch[:,20:65 ,25:50 ] = False
        print ('extracting flux H band XID 587 sinfoni')
    
    storage['Signal_mask'] = mask_catch
    
    
    if plot==1:
        plt.figure()
        plt.title('Selected Spaxels for 1D spectrum + Contours from 2D Gaus')
        plt.imshow(np.ma.array(data=storage['Median_stack_white'], mask=storage['Signal_mask'][0,:,:]), origin='low')
        plt.colorbar()

        shapes = storage['dim']
        x = np.linspace(0, shapes[1]-1, shapes[1])
        y = np.linspace(0, shapes[0]-1, shapes[0])
        x, y = np.meshgrid(x, y)

        data_fit = twoD_Gaussian((x,y), *storage['Median_stack_white_Center_data'])

        plt.contour(x, y, data_fit.reshape(shapes[0], shapes[1]), 8, colors='w')
   
    

    return storage



def astrometry_correction(storage):
    '''
    Correcting the Astrometry of the Cube. Fits a 2D Gaussian to the HST image and assumes that the 
    HST and Cube centroids are in the same location. 
    '''
    
    # Loading the HST image and the header
    img_file = '/Users/jansen/Google Drive/Astro/KMOS_SIN/HST_data/HST_'+storage['X-ray ID']+'.fits'
    
    img=fits.getdata(img_file)
    img_wcs= wcs.WCS(img_file).celestial
    hdr=fits.getheader(img_file)
    
    # Sie of the new image - same as the cube
    new_size = storage['Phys_size'] # size of the cube
    # Finding the pixel scale of the HST
    try:
        pixscale=abs(hdr['CD1_1']*3600)
    
    except:
        pixscale=abs(hdr['CDELT1']*3600)
        
    
    # Loading the Catalogue coordinates - Chris sometimes uses ra and sometimes RA
    Cat = storage['Cat_entry']
    try:
        Ra_opt = Cat['ra']
        Dec_opt = Cat['dec']
    
    except:
        Ra_opt = Cat['RA']
        Dec_opt = Cat['DEC']
    
    
    # Finding the position of the Galaxy in pix scale
    opt_world= np.array([[Ra_opt,Dec_opt]])
    opt_pixcrd = img_wcs.wcs_world2pix(opt_world, 0) # WCS transform
    opt_x= (opt_pixcrd[0,0]) # X pixel
    opt_y= (opt_pixcrd[0,1]) # Y pixel
    
    position = np.array([opt_x, opt_y])


    # Cutting out an image from the bigger image
    cutout = Cutout2D(img, position, new_size/pixscale, wcs=img_wcs,mode='partial')
       
    # Extracting the new image and the new wcs solution 
    img=(cutout.data).copy()
    img_wcs=cutout.wcs
    
    
    
    # To avoid weird things on side of the stamps
    img[np.isnan(img)] = 0
  
    # Finding the dimensions of the new image - need to create new XY grid for the 2D Gaussian 
    shapes = np.array(np.shape(img))
    
    # Finding the initial condition- location of the maximum of the image
    loc = np.where(img == img.max()) # Assumes that our AGN is the brightest thing in the image
    
    initial_guess = (np.max(img),loc[1][0],loc[0][0],5. , 5. ,0,0)
    
    # XID 522 in KMOS is just damn aweful 
    if storage['X-ray ID'] == 'XID_522':
        print( 'XID_522: Changing initial conditions for Astrometry corrections')
        initial_guess = (img[34,32],34,32,5. , 5. ,0,0)
        
        print (initial_guess)
        
        
    #print 'Initial guesses', initial_guess
    import scipy.optimize as opt
    
    x = np.linspace(0, shapes[1]-1, shapes[1])
    y = np.linspace(0, shapes[0]-1, shapes[0])
    x, y = np.meshgrid(x, y)   
    
    popt, pcov = opt.curve_fit(twoD_Gaussian, (x, y), img.ravel(), p0=initial_guess)
    er = np.sqrt(np.diag(pcov))
    print ('HST Cont', popt[1:3])
    print ('HST Cont er', er[1:3])
    popt_HST  = popt.copy()
    # Extrating the XY coordinates of the center. Will be plotted later to check accuracy
    center_C= popt[1:3]
    
    center_C1 = np.zeros(2)
    center_C1[0] = center_C[1]
    center_C1[1] = center_C[0]
    
    # Finding the Accurate HST based optical coordinates
    center_global =  img_wcs.wcs_pix2world(np.array([center_C]), 0)[0] 
    
    cube_center = storage['Median_stack_white_Center_data'][1:3]    

    # Old Header for plotting purposes
    Header_cube_old = storage['header'].copy()
    
    # New Header 
    Header_cube_new = storage['header'].copy()
    Header_cube_new['CRPIX1'] = cube_center[0]+1
    Header_cube_new['CRPIX2'] = cube_center[1]+1
    Header_cube_new['CRVAL1'] = center_global[0]
    Header_cube_new['CRVAL2'] = center_global[1]
    
    
    # Saving new coordinates and the new header
    storage['HST_center_glob'] = center_global
    storage['header'] = Header_cube_new
    
    ###############
    # DO check plots:    
    f  = plt.figure(figsize=(10,7))
    ax = f.add_subplot(121, projection=img_wcs)
    
    rms = np.std(img - img.mean(axis=0))
    
    ax.imshow((img), origin='low')#, vmin=0, vmax=3*rms)
    ax.set_autoscale_on(False)
    ax.plot(center_C[0], center_C[1], 'r*')
    
   
    cube_wcs= wcs.WCS(Header_cube_new).celestial  
    cont_map = storage['Median_stack_white']
    popt_cube = storage['Median_stack_white_Center_data']
    
    x = np.linspace(0, shapes[1]-1, shapes[1])
    y = np.linspace(0, shapes[0]-1, shapes[0])
    x, y = np.meshgrid(x, y)  
    
    
    data_fit = twoD_Gaussian((x,y), *popt_cube)
    
    ax.contour( data_fit.reshape(shapes[0], shapes[1]), levels=(max(data_fit)*0.68,max(data_fit)*0.98),transform= ax.get_transform(cube_wcs), colors='r')
    
    cube_wcs= wcs.WCS(Header_cube_old).celestial  
    cont_map = storage['Median_stack_white']
      
    ax.contour( data_fit.reshape(shapes[0], shapes[1]), levels=(max(data_fit)*0.68,max(data_fit)*0.98),transform= ax.get_transform(cube_wcs), colors='g')
      

    popt = cube_center
    cube_wcs= wcs.WCS(Header_cube_new).celestial
    
    ax = f.add_subplot(122, projection=cube_wcs)  
    ax.imshow(cont_map, vmin=0, vmax= cont_map[int(popt[1]), int(popt[0])], origin='low')
    
    ax.plot(popt[0], popt[1],'r*')
    
    ax.contour(img, transform= ax.get_transform(img_wcs), colors='red',levels=(popt_HST[0]*0.68,popt_HST[0]*0.95), alpha=0.9)
    
    ax.set_xlim(40,60)
    ax.set_ylim(40,57)
    
    return storage


def astrometry_correction_GAIA(storage):
    '''
    Correcting the Astrometry of the Cube. Fits a 2D Gaussian to the HST image and assumes that the 
    HST and Cube centroids are in the same location. 
    '''
    
    from astropy.table import Table, join, vstack, Column
    Gaia = Table.read(ph.MyPATH+'Catalogues/'+storage['X-ray ID']+'_Gaia.tbl' , format='ipac')
    
    
    
    
    
    cube_center = storage['Median_stack_white_Center_data'][1:3]    

    # Old Header for plotting purposes
    Header_cube_old = storage['header'].copy()
    
    # New Header 
    Header_cube_new = storage['header'].copy()
    Header_cube_new['CRPIX1'] = cube_center[0]+1
    Header_cube_new['CRPIX2'] = cube_center[1]+1
    Header_cube_new['CRVAL1'] = float(Gaia['ra'])
    Header_cube_new['CRVAL2'] = float(Gaia['dec'])
    
    center_global= np.array([float(Gaia['ra']), float(Gaia['dec'])])
    
    # Saving new coordinates and the new header
    storage['HST_center_glob'] = center_global
    storage['header'] = Header_cube_new
    
    
    
    return storage

#astrometry_correction(storage)

def D1_spectra_collapse(storage,z, band, plot, addsave=''):
    '''
    This function collapses the Cube to form a 1D spectrum of the galaxy
    '''
    # Loading the Signal mask - selects spaxel to collapse to 1D spectrum 
    mask_spax = storage['Signal_mask'].copy()
    # Loading mask of the sky lines an bad features in the spectrum 
    mask_sky_1D = storage['sky_clipped_1D'].copy()
    
    
    total_mask = mask_spax 
    # M
    flux = np.ma.array(data=storage['flux'].data, mask= total_mask) 
    
    D1_spectra = np.ma.sum(flux, axis=(1,2))
    D1_spectra = np.ma.array(data = D1_spectra.data, mask=mask_sky_1D)
    
    D1_s = np.ma.sum(np.ma.array(data=storage['flux'].data, mask= mask_spax) , axis=(1,2))
    
    wave= storage['obs_wave']
    
    if plot==1:
        plt.figure()
        plt.title('Collapsed 1D spectrum from D1_spectra_collapse fce')
        
        plt.plot(wave, D1_s, drawstyle='steps-mid', color='grey')
        plt.plot(wave, np.ma.array(data= D1_spectra, mask=storage['sky_clipped_1D']), drawstyle='steps-mid')
        '''
        OIIIa=  500.7*10/1e3*(1+z)
        OIIIb=  496.*10/1e3*(1+z)
        Hbeta=  485.*10/1e3*(1+z)
        Halpha = 6562.8/1e4*(1+z)
    
        y_line = np.array([0, max(D1_spectra)*1.2])
        if band=='OIII':
            plt.plot( np.array([OIIIa, OIIIa])  ,y_line, 'k--')
            plt.plot( np.array([OIIIb, OIIIb])  ,y_line, 'k--')
            plt.plot( np.array([Hbeta, Hbeta])  ,y_line, 'k--')
    
        if band =='H':
            plt.plot( np.array([Halpha, Halpha])  ,y_line, 'k--')
        '''
    
        plt.ylabel('Flux')
        plt.xlabel('Observed wavelength')
   
    
    storage['1D_spectrum'] = D1_spectra
    storage['1D_spectrum_er'] = STD_calc(wave/(1+z)*1e4,D1_spectra, band)* np.ones(len(D1_spectra))
    
    Save_spec = np.zeros((4,len(D1_spectra)))
    
    Save_spec[0,:] = wave
    Save_spec[1,:] = D1_spectra
    Save_spec[2,:] = storage['1D_spectrum_er'].copy()
    Save_spec[3,:] = mask_sky_1D
    
    ID = storage['X-ray ID']
    
    np.savetxt(ph.MyPATH+'Results_storage/Spectrums/'+ID+'_'+band+addsave+'.txt', Save_spec)
    
    return storage


def stack_sky(storage, Band,plot, spe_ma=np.array([], dtype=bool), expand=0):
    header = storage['header']
    ID = storage['X-ray ID']
    ######
    # Finding the center of the Object
    center =  storage['Median_stack_white_Center_data'][1:3].copy()
    
    ######
    # Defining other variable to be used
    wave = storage['obs_wave'].copy()#*1e4/(1+z) # Wavelength
    flux = storage['flux'].copy()   # Flux 
    
    shapes = storage['dim'] # The shapes of the image
    
    mask_nan = storage['flux'].mask.copy()
    
    mask_collapse = storage['flux'].mask.copy() # Mask to collapse 
    mask_collapse[:,:,:] = False
    
    
    
    header  = storage['header']
    try:
        arc = np.round(1.3/(header['CD2_2']*3600))
    
    except:
        arc = np.round(1.3/(header['CDELT2']*3600))
        
    
    for ix in range(shapes[0]):
        for iy in range(shapes[1]):
            dist = np.sqrt((ix- center[1])**2+ (iy- center[0])**2)
            if dist< arc:
                mask_collapse[:,ix,iy] = True
    
    try:
        arc = (header['CD2_2']*3600)
    except:
        arc = (header['CDELT2']*3600)
        
        
    if arc< 0.17:
        print ('Masking special corners')
        mask_collapse[:,:,0] = True
        mask_collapse[:,:,1] = True
        mask_collapse[:,:,2] = True
        mask_collapse[:,:,3] = True
        
        mask_collapse[:,:,-1] = True
        mask_collapse[:,:,-2] = True
        mask_collapse[:,:,-3] = True
        mask_collapse[:,:,-4] = True
        
        mask_collapse[:,0,:] = True
        mask_collapse[:,1,:] = True
        mask_collapse[:,2,:] = True
        mask_collapse[:,3,:] = True
        
        mask_collapse[:,-1,:] = True
        mask_collapse[:,-2,:] = True
        mask_collapse[:,-3,:] = True
        mask_collapse[:,-4,:] = True
    
    elif arc> 0.17:
        print ('Masking special corners')
        mask_collapse[:,:,0] = True
        mask_collapse[:,:,-1] = True  
        mask_collapse[:,0,:] = True       
        mask_collapse[:,-1,:] = True
        
    # For Hband sinfoni data I subtracted the frames from each other to imitate the sky frame. Therefore there are two negative objects in the cube. I am masking those as well.
    if (Band=='Hsin'):
        
        mask_collapse[:,76-10:76+8 ,83-10:83+10 ] = True #mask_collapse[:,2:10 ,0:12 ] = True #
        mask_collapse[:,21-10:21+8 ,27-10:27+10 ] = True
        print ('HB89 Hsin masking of negative objects')
        
    
    mask_collapse = np.logical_or(mask_collapse, mask_nan)
              
    storage['Sky_stack_mask'] = mask_collapse
    ######
    
    # Collapsing the spaxels into sky spectrum
    stacked_sky = np.ma.sum(np.ma.array(data=flux.data, mask=mask_collapse), axis=(1,2)) 
    
    # Creating a brand new mask to mask weird features in the sky_spectrum 
    std_mask = mask_collapse[:,7,7].copy()
    std_mask[:] = False
    
    weird = np.array([], dtype=int)
    
    if Band=='YJ':
        print ('Masking based on YJ band')
        # Actually masking the weird features: Edges + that 1.27 annoying bit
        weird = np.where((wave>1.25) &(wave<1.285))[0]
        weird = np.append(weird, np.where((wave<1.020))[0])
        weird = np.append(weird, np.where((wave>1.35))[0])
        weird = np.append(weird, np.where((wave>1.31168)&(wave< 1.3140375))[0])
    
        # Masking the stacked sky not to intefere with calculating std
        std_mask[weird] = True   
    
    elif  Band=='H':
        print ('Masking based on H band')
        weird = np.where((wave<1.45))[0]
        weird = np.append(weird, np.where((wave>1.85))[0])
        weird = np.append(weird, np.where((wave>1.78351) &(wave<1.7885))[0])
    
    elif Band=='K':
        print ('Masking based on K band')
        weird = np.where((wave<1.945))[0]
        weird = np.append(weird, np.where((wave>2.4))[0])
        
    elif Band=='Hsin':
        print ('Masking based on H band with Sinfoni')
        weird = np.where((wave>1.82))[0]
        weird = np.append(weird, np.where((wave<1.45))[0])
        print (len(wave), len(weird))
    
    elif Band=='Ysin':
        print ('Masking based on Y band with Sinfoni')
        weird = np.where((wave>1.35))[0]
        weird = np.append(weird, np.where((wave<1.11))[0])
        print (len(wave), len(weird))
    
    elif Band=='Ksin':
        print ('Masking based on K band with Sinfoni')
        weird = np.where((wave>2.4))[0]
        weird = np.append(weird, np.where((wave<1.945))[0])
        print (len(wave), len(weird))
        
    elif Band=='HKsin':
        print ('Masking based on HK band with Sinfoni')
        weird = np.where((wave>2.4))[0]
        weird = np.append(weird, np.where((wave<1.5))[0])
        print (len(wave), len(weird))
        
        
    
    if (ID=='2QZJ') & (Band=='Hsin'):
        print ('Masking based on 2QZJ Hsin')
        weird  = np.append(weird, [1301, 1302, 1303, 1304, 1305, 1306, 1307, 1308,1309,1310,1311,1312,1313,1368,1369,1370,1371])
        
    if (ID=='HB89') & (Band=='Hsin'):
        print ('Masking based on HB89 Hsin')
        weird  = np.append(weird, [1210, 1211, 1212, 1218, 1219, 1220])
        
    if (ID=='LBQS') & (Band=='Hsin'):
        print ('Masking based on LBQS Hsin')
        weird  = np.append(weird, [811,812,882,883,884,885,886,887,888,889, 1017,1018])
        
        
        
    
    weird = np.append(weird, spe_ma)
        
    std_mask[weird] = True

    stacked_sky_mask = np.ma.array(data = stacked_sky, mask=std_mask)
    
    
    # Masking sky lines 
    y = np.zeros_like(stacked_sky)
        
    clip_w = 1.5
    
    ssma = stacked_sky_mask.data[np.invert(stacked_sky_mask.mask)]
    print (len(ssma))
    low, hgh = conf(ssma)   
    
    
    sky = np.where((stacked_sky<y+ low*clip_w) | (stacked_sky> y + hgh*clip_w))[0]
    
    sky_clipped =  storage['flux'].mask.copy()
    sky_clipped = sky_clipped[:,7,7]
    sky_clipped[sky] = True          # Masking the sky features
    sky_clipped[weird] = True        # Masking the weird features
        
      
    # Storing the 1D sky line mask into cube to be used later
    mask_sky = storage['flux'].mask.copy()
    

    
    #if (storage['X-ray ID']=='HB89') & (Band=='Hsin'):
    #    sky_clipped[np.where((wave< 1.72344 ) & (wave > 1.714496))[0]] = False
    
    
        
    
    for ix in range(shapes[0]):
        for iy in range(shapes[1]):
            mask_sky[:,ix,iy] = sky_clipped
    
    # Storing the 1D and 3D sky masks
    storage['sky_clipped'] = mask_sky
    storage['sky_clipped_1D'] = sky_clipped 
    
    storage['Stacked_sky'] = stacked_sky
    
    
    #np.savetxt(ph.MyPATH+ID+'_Skyline_mask.txt', sky_clipped)
    

    
    if plot==1:
        mask_spax = storage['Signal_mask'].copy()
        mask_sky = storage['sky_clipped'].copy()
    
        total_mask = np.logical_or(mask_spax, mask_sky)
    
    
        flux_new = np.ma.array(data=storage['flux'].data, mask= total_mask) 
        D1_spectra_new = np.ma.sum(flux_new, axis=(1,2))
    
        flux_old = np.ma.array(data=storage['flux'].data, mask=mask_spax)
        D1_spectra_old = np.ma.sum(flux_old, axis=(1,2))
        #######
        # Plotting the spaxels used to assemble the sky
        
        plt.figure()
        plt.imshow(np.ma.array(data=storage['Median_stack_white'], mask=mask_collapse[0,:,:]), origin='low')
        plt.colorbar()
        plt.title('Removing the Galaxy')
    
    
        #######
        # 3 panel plot
        # 1 panel - The stacked spectrum
        # 2 panel - clipped sky lines
        # 3 panel - The collpsed 1D spectrum of the galaxy - No masking
        # 4 panel - The collpsed 1D spectrum of the galaxy - Final masking
        
        f, (ax1,ax2, ax3, ax4) = plt.subplots(4, sharex=True, figsize=(10,15))
        ax1.set_title('Stacked spectrum outside the galaxy')
    
        ax1.plot(wave, y+np.ones_like(stacked_sky)*hgh, 'g--')
        ax1.plot(wave, y+np.ones_like(stacked_sky)*low, 'g--')
        
        ax1.plot(wave, y+np.ones_like(stacked_sky)*hgh*clip_w, 'r--')
        ax1.plot(wave, y+np.ones_like(stacked_sky)*low*clip_w, 'r--')
        
        ax1.plot(wave, (stacked_sky), drawstyle='steps-mid', color='grey')
        ax1.plot(wave, np.ma.array(data=stacked_sky,mask=std_mask), drawstyle='steps-mid')
        ax1.set_ylabel('Sky spec')
        
        ax1.set_ylim(np.ma.min(np.ma.array(data=stacked_sky,mask=std_mask)), np.ma.max(np.ma.array(data=stacked_sky,mask=std_mask)))
    
        ax2.plot(wave, np.ma.array(data=stacked_sky, mask=sky_clipped), drawstyle='steps-mid')
        ax2.set_ylabel('Clipped sky')
        
        ax2.set_ylim(np.ma.min(np.ma.array(data=stacked_sky,mask=std_mask)), np.ma.max(np.ma.array(data=stacked_sky,mask=std_mask)))
    
        
        ax3.set_ylabel('Old spec')
        ax3.plot(wave, D1_spectra_old, drawstyle='steps-mid')
        
        ax4.plot(wave, D1_spectra_new, drawstyle='steps-mid')
        ax4.set_ylabel('New spec')
        
    
    return storage


def SNR_calc(flux, wave, solution, mode,z, mul=0):
    std = STD_calc(wave/(1+z)*1e4,flux, mode)
    wave = wave[np.invert(flux.mask)]
    flux = flux.data[np.invert(flux.mask)]
    
    wv_or = wave.copy()
    
    
    if (mode =='OIII') &(mul==0):
        center = solution.params['o3r_center'].value        
        fwhm = solution.params['o3r_fwhm'].value *2                
        use = np.where((wave< center+fwhm)&(wave> center-fwhm))[0]
        
    elif (mode =='OIII') & (mul==1):
        center = solution.params['o3rw_center'].value        
        fwhm = solution.params['o3rw_fwhm'].value *2                
        use = np.where((wave< center+fwhm)&(wave> center-fwhm))[0]
        
     
    elif (mode=='H') &(mul==0):
        center = solution.params['Ha_center'].value
        fwhm = solution.params['Ha_fwhm'].value *1.
        use = np.where((wave< center+fwhm)&(wave> center-fwhm))[0]
    
    elif (mode=='H') &(mul==1):
        center = solution.params['Han_center'].value
        fwhm = solution.params['Han_fwhm'].value *1.
        use = np.where((wave< center+fwhm)&(wave> center-fwhm))[0]
        
        #print center, fwhm, len(use), std
    
    elif (mode=='Hb') &(mul==0):
        center = solution.params['Hb_center'].value
        fwhm = solution.params['Hb_fwhm'].value *1
        use = np.where((wave< center+fwhm)&(wave> center-fwhm))[0]
    
    elif (mode=='Hb') &(mul==1):
        center = solution.params['Hbn_center'].value
        fwhm = solution.params['Hbn_fwhm'].value *1
        use = np.where((wave< center+fwhm)&(wave> center-fwhm))[0]
    
    elif (mode=='Hbs'):
        center = solution.params['Hb_center'].value
        fwhm = solution.params['Hb_fwhm'].value *1
        use = np.where((wave< center+fwhm)&(wave> center-fwhm))[0]
        
        #print center, fwhm, len(use), std
        
        
    flux_p = flux[use]
    wave = wave[use]
    
    import scipy.integrate as scpi
    try:
        Fit = scpi.simps( solution.eval(x=wave),wave )
        Dat = scpi.simps( flux_p,wave )
    
        rat = Fit/Dat
    
    except:
        rat = 1
    if rat<0:
        rat=0
    
    if (mul==1) & (mode=='H'):
        flux_l = flux_p - (solution.eval_components(x= wave)['linear'] + 
                           solution.eval_components(x= wave)['Haw_']+
                           solution.eval_components(x= wave)['Nr_']+
                           solution.eval_components(x= wave)['Nb_'] )
    if (mul==1) & (mode=='OIII'):
        flux_l = flux_p - (solution.eval_components(x= wave)['linear'] + 
                           solution.eval_components(x= wave)['Hbw_']+
                           solution.eval_components(x= wave)['o3bn_']+
                           solution.eval_components(x= wave)['o3bw_']
                            )
        
        
        
        
    else:
        flux_l = flux_p -  solution.eval_components(x= wave)['linear']
        
        
    
    n = len(use)
    if std==0:
        SNR=0
    
    
    else:
        SNR = (sum(flux_l)/np.sqrt(n)) * (1./std)
        #print SNR
       
        if SNR < 0:
            SNR=0
    #print 'std ',std, ', #points ', n,  'SNR ', SNR, 'sum flux ', sum(flux_l)
    return SNR, rat
    
            

def fitting_collapse_OIII(storage,z, plot):
    wave = storage['obs_wave'].copy()
    flux = storage['1D_spectrum'].copy()
    error = storage['1D_spectrum_er'].copy()
    
    ID = storage['X-ray ID']
    
    
    fit_loc = np.where((wave>4400*1e4/(1+z))&(wave<5300*1e4/(1+z)))[0]
    #if ID=='HB89':
    #    out = emfit.fitting_OIIIbkp_Hbeta_qso_mul(wave,flux,error,z)
    #else:
    
    if ID=='LBQS':
        out = emfit.fitting_OIII_Hbeta_qso_mul(wave,flux, error,z, hbw=4870.)
    else:
        out = emfit.fitting_OIII_Hbeta_qso_mul(wave,flux, error,z, hbw=4861.)
    
    
    
    outs, chi2s = emfit.fitting_OIII_sig(wave,flux,error,z)
    
    SNR,dat = SNR_calc(flux, wave, outs, 'OIII',z)
    
    
    print( 'SNR of the line is ', SNR, ' and recovers ', dat, ' of the flux')
    #print( 'The total flux of the OIII is ', flux_measure_ind(out,wave, 'OIII', use='tot' ))
    
    if plot==1:
        f = plt.figure()
        ax1 = f.add_axes([0.125, 0.2, 0.8,0.7])

        axres = f.add_axes([0.125, 0.1, 0.8, 0.1])
        
        
    
        ax1.set_title('Collapsed 1D spectrum')
        ax1.set_xlabel('Rest Wavelegth (ang)')
        ax1.set_ylabel('Flux')
    
        from Tools_plotting import plotting_OIII_Hbeta
        plotting_OIII_Hbeta(wave, flux, ax1, out, 'mul',z)
        ax1.plot(wave[fit_loc], outs.eval(x=wave[fit_loc]), 'k--')
        
        #axres.plot(wave[fit_loc], out.eval(x=wave[fit_loc]), drawstyle='steps-mid' )
        #axres.set_xlim(4700, 5100)
    
    storage['OIII_fit_collapse']= out   
    storage['1D_fit_OIII_sig'] = outs
    
    storage['z_fit'] = (out.params['o3rn_center']*1e4/5008.)-1
    
    print ('Old redshift ', z, ' new z ', (out.params['o3rn_center']*1e4/5008.)-1 )

    
    storage['1D_fit_OIII_mul'] = out
    storage['1D_fit_OIII_sig'] = outs
    storage['1D_fit_OIII_SNR'] = SNR
    
    return storage


def fitting_collapse_Halpha(storage,z, plot, broad = 1, cont=1):
    wave = storage['obs_wave'].copy()
    flux = storage['1D_spectrum'].copy()
    error = storage['1D_spectrum_er'].copy()
    
    ID = storage['X-ray ID']
    
    fl = flux.data
    ms = flux.mask
    
    #SII_ms = ms.copy()
    #SII_ms[:] = False
    #SII_ms[np.where(((wave*1e4/(1+z))<6741)&((wave*1e4/(1+z))> 6712))[0]] = True
    
    #msk = np.logical_or(SII_ms,ms)  
    msk = ms
    
    flux = np.ma.array(data=fl, mask = msk)
    
    fit_loc = np.where(((wave*1e4/(1+z))>6562.8-100)&((wave*1e4/(1+z))<6562.8+100))[0]
    
    out, chi2 = emfit.fitting_Halpha_mul_testing(wave,flux,error,z, broad=broad, cont=cont)
    
    outs = emfit.fitting_Halpha_sig(wave,flux,error,z)
    
    SNR,chi2 = SNR_calc(flux, wave, outs, 'H',z)
    
    SNR,chi2 = SNR_calc(flux, wave, out, 'H',z, mul=1)
    
    
    if plot==1:    
        f, (ax1) = plt.subplots(1)
    
        ax1.set_title('Collapsed 1D spectrum')
        ax1.set_xlabel('Rest Wavelegth (ang)')
        ax1.set_ylabel('Flux')
        
        from Tools_plotting import plotting_Halpha
    
        plotting_Halpha(wave, flux, ax1, out, 'mul',z)
        ax1.plot(wave[fit_loc]/(1+z)*1e4, outs.eval(x=wave[fit_loc]), 'k--')
    
    
    #Width_narrow = out.params['Han_fwhm'].value*1e4/(1+z)/6562.8*2.9979e5
    #Width_narrow_Er = out.params['Han_fwhm'].stderr/6562.8*2.9979e5
    
    print ('The width of a single Gauss', outs.params['Ha_fwhm'].value*1e4/(1+z)/6562.8*2.9979e5)
    
    #print ('The width of a mul Gauss', out.params['Haw_fwhm'].value*1e4/(1+z)/6562.8*2.9979e5, out.params['Han_fwhm'].value*1e4/(1+z)/6562.8*2.9979e5)
    
    print (' ')
    print ('Flux Halpha narrow ', flux_measure_ind(out,wave, 'H', use='BPT'))
    print ('Flux Halpha broad ', flux_measure_ind(out,wave, 'H', use='broad'))
    print ('Flux NII ', flux_measure_ind(out,wave, 'NII', use='tot'))
    print (' ')
    
    storage['Halpha_fit_collapse']= out
      
    storage['1D_fit_Halpha_mul'] = out
    storage['1D_fit_Halpha_sig'] = outs
    storage['1D_fit_Halpha_SNR'] = SNR    
    
    return storage





def Spaxel_fit_wrap_sig(storage, Line_info, obs_wv, flx_spax_m, error, mode,i,j ,broad):
    obs_wave = storage['obs_wave']
    
    z = storage['z_guess']
    
    if mode=='H':
        out_list = []
        chi_list = np.array([])
        Hal_cm = 6562.8*(1+z)/1e4
        
        suc=1
        
        D_out = storage['1D_fit_Halpha_mul']
        loc = D_out.params['Han_center'].value 
            
              
        dec = D_out
            
        
        try:
            out,chi2 = emfit.fitting_Halpha_mul_testing(obs_wv,flx_spax_m ,error,z,broad=broad, decompose=  dec, init_sig=250.,  wvnet=loc)
            out_list.append(out)
            chi_list = np.append(chi_list, chi2)
            
            out,chi2 = emfit.fitting_Halpha_mul_testing(obs_wv,flx_spax_m ,error,z,broad=broad, decompose=  dec,init_sig=370.,   wvnet=loc)
            out_list.append(out)
            chi_list = np.append(chi_list, chi2)
            
            out,chi2 = emfit.fitting_Halpha_mul_testing(obs_wv,flx_spax_m ,error,z,broad=broad, decompose=  dec, init_sig=450.,  wvnet=loc)
            out_list.append(out)
            chi_list = np.append(chi_list, chi2)
            
            
            best = np.argmin(chi_list)
            out = out_list[best]
            
            
            SNR,dat = SNR_calc(flx_spax_m,obs_wv, out, 'H',z, mul=1)
            
            Line_info[3,i,j] = SNR  
            Line_info[4,i,j] = dat
            
                    
                    
            if SNR>3:
                Line_info[0,i,j] =   flux_measure_ind(out,obs_wave, 'H', use='BPT')  #out.params['Han_height'].value
                Line_info[1,i,j] = -((loc-out.params['Han_center'].value)/Hal_cm)*2.9979e5
                Line_info[2,i,j] = (out.params['Han_fwhm'].value/Hal_cm)*2.9979e5  
                Line_info[5,i,j] =  out.params['Haw_amplitude'].value
                
        
            else:
                Line_info[0,i,j] = np.nan
                Line_info[1,i,j] = np.nan
                Line_info[2,i,j] = np.nan
                Line_info[5,i,j] = out.params['Haw_amplitude'].value
        
        
        except:
            print ('Spaxel fit fail')
            suc=0                    
            Line_info[0,i,j] = np.nan
            Line_info[1,i,j] = np.nan
            Line_info[2,i,j] = np.nan
            Line_info[3,i,j] = np.nan
            SNR=0
            out=1
                               
    elif mode =='OIII':
        out_old = storage['1D_fit_OIII_sig']
        
        
        OIII_cm = 5006.9*(1+z)/1e4
        width = (out_old.params['o3r_fwhm'].value/OIII_cm)*2.9979e5
        
        out_list = []
        chi_list = np.array([])
        
        try:
            suc = 1
            '''
            out, chi2 = emfit.fitting_OIII_sig(obs_wv,flx_spax_m ,error,z, init_sig=width, init_offset=0.)
            out_list.append(out)
            chi_list = np.append(chi_list, chi2)
            
            out, chi2 = emfit.fitting_OIII_sig(obs_wv,flx_spax_m ,error,z, init_sig=width, init_offset=-2.5/1e4*(1+z))
            out_list.append(out)
            chi_list = np.append(chi_list, chi2)
            
            out, chi2 = emfit.fitting_OIII_sig(obs_wv,flx_spax_m ,error,z, init_sig=width, init_offset=2.5/1e4*(1+z))
            out_list.append(out)
            chi_list = np.append(chi_list, chi2)
            '''
            
            out, chi2 = emfit.fitting_OIII_sig(obs_wv,flx_spax_m ,error,z, init_sig=width*1.2)
            out_list.append(out)
            chi_list = np.append(chi_list, chi2)
            
            out, chi2 = emfit.fitting_OIII_sig(obs_wv,flx_spax_m ,error,z, init_sig=width*1.5)
            out_list.append(out)
            chi_list = np.append(chi_list, chi2)
            
            out, chi2 = emfit.fitting_OIII_sig(obs_wv,flx_spax_m ,error,z, init_sig=width*0.66)
            out_list.append(out)
            chi_list = np.append(chi_list, chi2)
            
            out, chi2 = emfit.fitting_OIII_sig(obs_wv,flx_spax_m ,error,z, init_sig=width*0.4)
            out_list.append(out)
            chi_list = np.append(chi_list, chi2)
            
            out, chi2 = emfit.fitting_OIII_sig(obs_wv,flx_spax_m ,error,z, init_sig=width)
            out_list.append(out)
            chi_list = np.append(chi_list, chi2)
            
            
            best = np.argmin(chi_list)
            out = out_list[best]
                    
            D_out = storage['1D_fit_OIII_sig']
            loc = D_out.params['o3r_center'].value 
            
            SNR,dat = SNR_calc(flx_spax_m,obs_wv, out, 'OIII',z)
            Line_info[3,i,j] = SNR  
            Line_info[4,i,j] = dat
                    
                    
            if SNR>3:
                Line_info[0,i,j] = flux_measure_ind(out,obs_wave, 'OIIIs', use='tot')
                Line_info[1,i,j] = -((loc- out.params['o3r_center'].value)/OIII_cm)*2.9979e5
                Line_info[2,i,j] = (out.params['o3r_fwhm'].value/OIII_cm)*2.9979e5
                        
                        
            else:
                Line_info[0,i,j] = np.nan
                Line_info[1,i,j] = np.nan
                Line_info[2,i,j] = np.nan
                   
        except:
            print ('Spaxel fit fail')
            suc=0                    
            Line_info[0,i,j] = np.nan
            Line_info[1,i,j] = np.nan
            Line_info[2,i,j] = np.nan
            Line_info[3,i,j] = np.nan
            SNR = 0
            out =1 
    
    return SNR, Line_info,out,suc


def Spaxel_fit_sig(storage, mode, plot, sp_binning, localised = 0, broad=1, instrument='KMOS', add=''):
    flux = storage['flux'].copy()
    Mask= storage['sky_clipped_1D']
    shapes = storage['dim']
    
    ThD_mask = storage['sky_clipped'].copy()
    z = storage['z_guess']
    wv_obs = storage['obs_wave'].copy()
    
    Residual = np.zeros_like(flux).data
    Model = np.zeros_like(flux).data
        
      
    ms = Mask.copy()
    Spax_mask = storage['Sky_stack_mask'][0,:,:]
    
    if (storage['X-ray ID']=='XID_587'):
        print ('Masking special corners')
        Spax_mask[:,0] = False
        Spax_mask[:,1] = False
        
        Spax_mask[:,-1] = False
        Spax_mask[:,-2] = False
        
        Spax_mask[0,:] = False
        Spax_mask[1,:] = False
        
        Spax_mask[-1,:] = False
        Spax_mask[-2,:] = False
        
    
    if mode =='H':
        SII_ms = ms.copy()
        SII_ms[:] = False
        SII_ms[np.where((wv_obs<6741.*(1+z)/1e4)&(wv_obs> 6712*(1+z)/1e4))[0]] = True
    
        msk = np.logical_or(SII_ms, ms)    
    
    elif mode=='OIII':
        msk = ms  
    
    ID = storage['X-ray ID']    
    Line_info = np.zeros((6, shapes[0], shapes[1]))
    
    if plot==1:
        f, (ax1) = plt.subplots(1)
        
        ax1.set_xlabel('Rest Wavelegth (ang)')
        ax1.set_ylabel('Flux')
        
        if sp_binning=='Individual':
            
            Spax = PdfPages(ph.MyPATH+'Graphs/Spax_fit/Spaxel_'+ID+'_'+mode+add+'.pdf')
            
            
        
        elif sp_binning=='Nearest':
            Spax = PdfPages(ph.MyPATH+'Graphs/Spax_fit/Spaxel_Nearest_'+ID+'_'+mode+add+'.pdf')
            
            
    
    
    import Tools_plotting as emplot
    #############################
    # Binning Spaxel Fitting
    #############################
    if sp_binning== 'Nearest':
        header = storage['header']
        
        try:
            arc = (header['CD2_2']*3600)
        
        except:
            arc = (header['CDELT2']*3600)
        
        Line_info[:,:,:] = np.nan
        
        if arc> 0.17:            
            upper_lim = 2            
            step = 1
            
            if localised==1:
                popt =  storage['Median_stack_white_Center_data'][1:3].copy()
                
                x = np.linspace(popt[0]-4, popt[0]+4, 9) 
                x =np.array(x, dtype=int)
                y = np.linspace(popt[1]-4, popt[1]+4, 9) 
                y =np.array(y, dtype=int)
            
            else:
                x = range(shapes[0]-upper_lim)
                y = range(shapes[1]-upper_lim)               
            
        elif arc< 0.17:            
            upper_lim = 3 
            step = 2
            
            if localised==1:
                popt =  storage['Median_stack_white_Center_data'][1:3].copy()
                
                x = np.linspace(popt[1]-12, popt[1]+12, 25) 
                x =np.array(x, dtype=int)
                y = np.linspace(popt[0]-12, popt[0]+12, 25) 
                y =np.array(y, dtype=int)
            
                               
            else:
                x = range(shapes[0]-upper_lim)
                y = range(shapes[1]-upper_lim)
        
        
        #x = np.array([20])
        #y = np.array([20])
        for i in x: #progressbar.progressbar(x):
            i= i+step
            print (i,'/',len(x))
            for j in y:
                
                j=j+step
                #print i,j
                Spax_mask_pick = ThD_mask.copy()
                Spax_mask_pick[:,:,:] = True
                Spax_mask_pick[:, i-step:i+upper_lim, j-step:j+upper_lim] = False
                
                #Spax_mask_pick= np.logical_or(Spax_mask_pick, ThD_mask)
                flx_spax_t = np.ma.array(data=flux,mask=Spax_mask_pick)
                
                flx_spax = np.ma.median(flx_spax_t, axis=(1,2))                
                flx_spax_m = np.ma.array(data = flx_spax.data, mask=msk)                
                error = STD_calc(wv_obs/(1+z)*1e4,flx_spax, mode)* np.ones(len(flx_spax))
          
                SNR, Line_info,out, suc = Spaxel_fit_wrap_sig(storage, Line_info, wv_obs, flx_spax_m, error, mode,i,j, broad )
                
                if out !=1:
                    try:
                        
                        out = out[0]
                    except:
                        out=out
                
                
                    Residual[:,i,j] = flx_spax.data - out.eval(x=wv_obs)
                    Model[:,i,j] =  out.eval(x=wv_obs)   
                
                               
                prhdr = storage['header']
                hdu = fits.PrimaryHDU(Line_info, header=prhdr)
                hdulist = fits.HDUList([hdu])
    
                if mode == 'OIII':  
                    hdulist.writeto(ph.MyPATH+'Results_storage/OIII/'+ID+'_'+sp_binning+'_spaxel_fit'+add+'.fits', overwrite=True)
                    
                    
        
        
                elif mode == 'H':
                    hdulist.writeto(ph.MyPATH+'Results_storage/Halpha/'+ID+'_'+sp_binning+'_spaxel_fit'+add+'.fits', overwrite=True)
                    
                    
                        
                
                if (plot==1) &(suc==1):
                    
                    if (Spax_mask[i,j]==True) :
                        ax1.set_title('Spaxel '+str(i)+', '+str(j)+' in Obj with SNR = "%.3f' % SNR )
                        
                                              
                        ax1.plot(wv_obs, flx_spax_m.data, color='grey', drawstyle='steps-mid')                       
                        ax1.plot(wv_obs[np.invert(flx_spax_m.mask)], flx_spax_m.data[np.invert(flx_spax_m.mask)], drawstyle='steps-mid')                   
                        ax1.plot(wv_obs, out.eval(x=wv_obs), 'r--')
        
                        if mode=='H':
                            ax1.set_xlim(6400.*(1+z)/1e4, 6700.*(1+z)/1e4)
                            
                            ax1.plot(wv_obs, out.eval_components(x=wv_obs)['Han_'], color='orange', linestyle='dashed')
                            ax1.plot(wv_obs, out.eval_components(x=wv_obs)['Haw_'], color='blue', linestyle='dashed')
                            ax1.plot(wv_obs, out.eval_components(x=wv_obs)['Nr_'], color='green', linestyle='dashed')
                            ax1.plot(wv_obs, out.eval_components(x=wv_obs)['Nb_'], color='limegreen', linestyle='dashed')
                    
                            Hal_cm = 6562.*(1+z)/1e4
                    
                            #print (out.params['Han_fwhm'].value/Hal_cm)*2.9979e5 
                    
                        elif mode=='OIII':
                            cen = out.params['o3r_center'].value
                            wid = out.params['o3r_fwhm'].value
                            use = np.where((wv_obs< cen+wid)&(wv_obs> cen-wid))[0]
                            
                            ax1.plot(wv_obs, out.eval(x=wv_obs), 'k--')
                            ax1.plot(wv_obs[use], out.eval_components(x= wv_obs[use])['o3r_'], 'r--')
                            ax1.set_xlim(4900.*(1+z)/1e4, 5100.*(1+z)/1e4)
                        
                        #Spax.savefig()    
                        ax1.clear()
    #############################
    # Individual Spaxel Fitting
    #############################
    elif sp_binning== 'Individual':
        
        header = storage['header']
        
        try:
            arc = (header['CD2_2']*3600)
        
        except:
            arc = (header['CDELT2']*3600)
        
        Line_info[:,:,:] = np.nan
        
        if arc> 0.17:            
            upper_lim = 2            
            step = 1
            
            if localised==1:
                popt =  storage['Median_stack_white_Center_data'][1:3].copy()
                
                x = np.linspace(popt[0]-4, popt[0]+4, 9) 
                x =np.array(x, dtype=int)
                y = np.linspace(popt[1]-4, popt[1]+4, 9) 
                y =np.array(y, dtype=int)
            
            else:
                x = range(shapes[0]-upper_lim)
                y = range(shapes[1]-upper_lim)               
            
        elif arc< 0.17:            
            upper_lim = 3 
            step = 2
            
            if localised==1:
                popt =  storage['Median_stack_white_Center_data'][1:3].copy()
                
                x = np.linspace(popt[1]-12, popt[1]+12, 25) 
                x =np.array(x, dtype=int)
                y = np.linspace(popt[0]-12, popt[0]+12, 25) 
                y =np.array(y, dtype=int)
            
                #x = np.linspace(35, 65, 31) 
                #x =np.array(x, dtype=int)
                #y = np.linspace(35, 65, 31) 
                #y =np.array(y, dtype=int)
                #x = np.array([45])-2
                #y = np.array([47])-2
                
            else:
                x = range(shapes[0]-upper_lim)
                y = range(shapes[1]-upper_lim)
        
        
        for i in x: #progressbar.progressbar(x):
            i= i+step
            print (i,'/',len(x))
            for j in y:
                
                
                flx_spax = flux[:,i,j]
                flx_spax_m = np.ma.array(data=flx_spax, mask= msk)
                error =   STD_calc(wv_obs*(1+z)/1e4,flx_spax, mode)* np.ones(len(flx_spax))
            
                SNR, Line_info,out, suc = Spaxel_fit_wrap_sig(storage, Line_info, wv_obs, flx_spax_m, error, mode,i,j, broad )
                
                prhdr = storage['header']
                hdu = fits.PrimaryHDU(Line_info, header=prhdr)
                
                hdulist = fits.HDUList([hdu])
    
                if mode == 'OIII':            
                    hdulist.writeto(ph.MyPATH+'Results_storage/OIII/'+ID+'_'+sp_binning+'_spaxel_fit'+add+'.fits', overwrite=True)
        
        
                elif mode == 'H':     
                    hdulist.writeto(ph.MyPATH+'Results_storage/Halpha/'+ID+'_'+sp_binning+'_spaxel_fit'+add+'.fits', overwrite=True)
        
                if (plot==1) &(suc==1): 
                    if (Spax_mask[i,j]==True) :
                        ax1.set_title('Spaxel '+str(i)+', '+str(j)+' in Obj with SNR = "%.3f' % SNR )
                        
                                              
                        ax1.plot(wv_obs, flx_spax_m.data, color='grey', drawstyle='steps-mid')                       
                        ax1.plot(wv_obs[np.invert(flx_spax_m.mask)], flx_spax_m.data[np.invert(flx_spax_m.mask)], drawstyle='steps-mid')                   
                        ax1.plot(wv_obs, out.eval(x=wv_obs), 'r--')
        
                        if mode=='H':
                            ax1.set_xlim(6400.*(1+z)/1e4, 6700.*(1+z)/1e4)
                            
                            ax1.plot(wv_obs, out.eval_components(x=wv_obs)['Han_'], color='orange', linestyle='dashed')
                            ax1.plot(wv_obs, out.eval_components(x=wv_obs)['Haw_'], color='blue', linestyle='dashed')
                            ax1.plot(wv_obs, out.eval_components(x=wv_obs)['Nr_'], color='green', linestyle='dashed')
                            ax1.plot(wv_obs, out.eval_components(x=wv_obs)['Nb_'], color='limegreen', linestyle='dashed')
                    
                            Hal_cm = 6562.*(1+z)/1e4
                    
                            #print (out.params['Han_fwhm'].value/Hal_cm)*2.9979e5 
                    
                        elif mode=='OIII':
                            cen = out.params['o3r_center'].value
                            wid = out.params['o3r_fwhm'].value
                            use = np.where((wv_obs< cen+wid)&(wv_obs> cen-wid))[0]
                            
                            ax1.plot(wv_obs, out.eval(x=wv_obs), 'k--')
                            ax1.plot(wv_obs[use], out.eval_components(x= wv_obs[use])['o3r_'], 'r--')
                            ax1.set_xlim(4900.*(1+z)/1e4, 5100.*(1+z)/1e4)
                        
                        Spax.savefig()    
                        ax1.clear()
    if plot==1:
        Spax.close()
    
    prhdr = storage['header']
    hdu = fits.PrimaryHDU(Line_info, header=prhdr)
    hdulist = fits.HDUList([hdu])
    
    hdu_res = fits.PrimaryHDU(Residual, header=prhdr)
    hdulist_res = fits.HDUList([hdu_res])
    
    hdu_mod = fits.PrimaryHDU(Model, header=prhdr)
    hdulist_mod = fits.HDUList([hdu_mod])
    
    if mode == 'OIII':
        hdulist.writeto(ph.MyPATH+'Results_storage/OIII/'+ID+'_'+sp_binning+'_spaxel_fit'+add+'.fits', overwrite=True)
        
    elif mode == 'H':
        hdulist.writeto(ph.MyPATH+'Results_storage/Halpha/'+ID+'_'+sp_binning+'_spaxel_fit'+add+'.fits', overwrite=True)
            
  
        
        hdulist_res.writeto(ph.MyPATH+'Results_storage/Halpha/'+ID+'_'+sp_binning+'_spaxel_fit'+add+'_res.fits', overwrite=True)
        hdulist_mod.writeto(ph.MyPATH+'Results_storage/Halpha/'+ID+'_'+sp_binning+'_spaxel_fit'+add+'_mod.fits', overwrite=True)
        
    return storage

import Fitting_tools as emfit

def Spaxel_fit_wrap_mul(storage, Line_info, obs_wv, flx_spax_m, error, Residual, mode,i,j ,broad):
    
    z = storage['z_fit']
                               
    if mode =='OIII':
        outo = storage['1D_fit_OIII_mul']
        
        
        OIII_cm = 5006.9*(1+z)/1e4
        #width = (outo.params['o3rn_fwhm'].value/OIII_cm)*2.9979e5
        
        out_list = []
        chi_list = np.array([])
        
        if 1==1:
            suc = 1
            
            out, chi2 = emfit.fitting_OIII_Hbeta_qso_mul(obs_wv,flx_spax_m, error,z, decompose=outo, chir=1)
            
            out_list.append(out)
            chi_list = np.append(chi_list, chi2)
            
            
            best = np.argmin(chi_list)
            out = out_list[0]
                    
            D_out = storage['1D_fit_OIII_mul']
            loc = D_out.params['o3rn_center'].value 
            
            SNR,dat = SNR_calc(flx_spax_m,obs_wv, out, 'OIII',z, mul=1)
            Line_info[9,i,j] = SNR  
                    
                    
            if SNR>3:
                Line_info[0,i,j] = out.params['o3rn_center'].value
                Line_info[1,i,j] = out.params['o3rn_sigma'].value
                Line_info[2,i,j] = out.params['o3rn_amplitude'].value
                
                Line_info[3,i,j] = out.params['o3rw_center'].value
                Line_info[4,i,j] = out.params['o3rw_sigma'].value
                Line_info[5,i,j] = out.params['o3rw_amplitude'].value
                
                Line_info[6,i,j] = out.params['Hbn_center'].value
                Line_info[7,i,j] = out.params['Hbn_sigma'].value
                Line_info[8,i,j] = out.params['Hbn_amplitude'].value
                
                Line_info[8,i,j] = out.params['Hbw_amplitude'].value
                Line_info[9,i,j] = out.params['Hbw_center'].value
                Line_info[10,i,j] = out.params['Hbw_sigma'].value
                Line_info[11,i,j] = out.params['Hbw_a1'].value
                Line_info[12,i,j] = out.params['Hbw_a2'].value   

                Line_info[13,i,j] = out.params['slope'].value
                Line_info[14,i,j] = out.params['intercept'].value
                
                Line_info[15,i,j] = SNR
                
                Residual[:,i,j] = flx_spax_m.data - out.eval(x=obs_wv)
                
                                    
            else:
                Line_info[:14,i,j] = np.nan
                
                   
        else:
            print ('Spaxel fit fail')
            suc=0                    
            Line_info[:,i,j] = np.nan
            
            SNR = 0
            out =1 
    
    return SNR, Line_info,out,suc, Residual


def Spaxel_fit_mul(storage, mode, plot, sp_binning, localised = 0, broad=1, instrument='KMOS', add=''):
    flux = storage['flux'].copy()
    Mask= storage['sky_clipped_1D']
    shapes = storage['dim']
    
    ThD_mask = storage['sky_clipped'].copy()
    z = storage['z_guess']
    wv_obs = storage['obs_wave'].copy()
        
      
    ms = Mask.copy()
    Spax_mask = storage['Sky_stack_mask'][0,:,:]  
    
    if mode=='OIII':
        msk = ms  
    
    ID = storage['X-ray ID']    
    Line_info = np.zeros((16, shapes[0], shapes[1]))
    
    Residual = flux.data.copy()
    
    Residual_OIII = Residual.copy()
    
    if plot==1:
        f, (ax1) = plt.subplots(1)
        
        ax1.set_xlabel('Rest Wavelegth (ang)')
        ax1.set_ylabel('Flux')
            
        
        if sp_binning=='Nearest':
            Spax = PdfPages(ph.MyPATH+'Graphs/Spax_fit/Spaxel_Nearest_'+ID+'_'+mode+add+'.pdf')
            
            
        elif sp_binning=='Individual':
            Spax = PdfPages(ph.MyPATH+'Graphs/Spax_fit/Spaxel_Nearest_'+ID+'_'+mode+add+'.pdf')
            
    
    
    
    import Tools_plotting as emplot
    #############################
    # Binning Spaxel Fitting
    #############################
    if 1==1:
        header = storage['header']
        
        try:
            arc = (header['CD2_2']*3600)
        
        except:
            arc = (header['CDELT2']*3600)
        
        Line_info[:,:,:] = np.nan
        
        if arc> 0.17:            
            upper_lim = 2            
            step = 1
            
            if localised==1:
                popt =  storage['Median_stack_white_Center_data'][1:3].copy()
                
                x = np.linspace(popt[0]-4, popt[0]+4, 9) 
                x =np.array(x, dtype=int)
                y = np.linspace(popt[1]-4, popt[1]+4, 9) 
                y =np.array(y, dtype=int)
            
            else:
                x = range(shapes[0]-upper_lim)
                y = range(shapes[1]-upper_lim)               
            
        elif arc< 0.17:            
            upper_lim = 3 
            step = 2
            
            if localised==1:
                popt =  storage['Median_stack_white_Center_data'][1:3].copy()
                
                x = np.linspace(popt[1]-12, popt[1]+12, 25) 
                x =np.array(x, dtype=int)
                y = np.linspace(popt[0]-12, popt[0]+12, 25) 
                y =np.array(y, dtype=int)
                

                
            else:
                x = range(shapes[0]-upper_lim)
                y = range(shapes[1]-upper_lim)
        
        #import progressbar
        #x = np.array([47])
        #y = np.array([50])
        for i in x: #progressbar.progressbar(x):
            i= i+step
            print (i-x[0]+1-step,'/',len(x))
            for j in y:
                
                j=j+step
                if sp_binning=='Nearest':
                    
                    #print i,j
                    Spax_mask_pick = ThD_mask.copy()
                    Spax_mask_pick[:,:,:] = True
                    Spax_mask_pick[:, i-step:i+upper_lim, j-step:j+upper_lim] = False
                    
                    #Spax_mask_pick= np.logical_or(Spax_mask_pick, ThD_mask)
                    flx_spax_t = np.ma.array(data=flux,mask=Spax_mask_pick)              
                    flx_spax = np.ma.median(flx_spax_t, axis=(1,2)) 
                
                elif sp_binning =='Individual':
                    flx_spax = flux[:,i,j].copy()
                    
                    
                flx_spax_m = np.ma.array(data = flx_spax.data, mask=msk)                
                error = STD_calc(wv_obs/(1+z)*1e4,flx_spax, mode)* np.ones(len(flx_spax))
                
                if broad==1:
                    broad= storage['1D_fit_OIII_mul']
          
                SNR, Line_info,out, suc, Residual = Spaxel_fit_wrap_mul(storage, Line_info, wv_obs, flx_spax_m, error, Residual , mode,i,j, broad )
                
                Residual_OIII[:,i,j] = flx_spax_m.data - (out.eval_components(x=wv_obs)['linear'] + out.eval_components(x=wv_obs)['Hbn_'] +  out.eval_components(x=wv_obs)['Hbw_']    ) 
                
                prhdr = storage['header']
                hdu = fits.PrimaryHDU(Line_info, header=prhdr)
                hdulist = fits.HDUList([hdu])
                
                prhdrr = storage['header']
                hdur = fits.PrimaryHDU(Residual, header=prhdrr)
                hdulistr = fits.HDUList([hdur])
                
                prhdrro = storage['header']
                hduro = fits.PrimaryHDU(Residual_OIII, header=prhdrr)
                hdulistro = fits.HDUList([hduro])
    
                if mode == 'OIII':
                    hdulist.writeto(ph.MyPATH+'Results_storage/OIII/'+ID+'_'+sp_binning+'_spaxel_fit.fits', overwrite=True)
                    hdulistr.writeto(ph.MyPATH+'Results_storage/OIII/'+ID+'_'+sp_binning+'_spaxel_residual.fits', overwrite=True)
                    hdulistro.writeto(ph.MyPATH+'Results_storage/OIII/'+ID+'_'+sp_binning+'_spaxel_residual_OIII.fits', overwrite=True)
                
                if (plot==1) &(suc==1):
                    
                    if (Spax_mask[i,j]==True) :
                        ax1.set_title('Spaxel '+str(i)+', '+str(j)+' in Obj with SNR = "%.3f' % SNR )
                                                                    
                        emplot.plotting_OIII_Hbeta(wv_obs, flx_spax_m, ax1, out, 'mul',z, title=0)
     
             
                        Spax.savefig()    
                        ax1.clear()
    
    if plot==1:
        Spax.close()
    
    prhdr = storage['header']
    hdu = fits.PrimaryHDU(Line_info, header=prhdr)
    hdulist = fits.HDUList([hdu])
    
    if mode == 'OIII':
        hdulist.writeto(ph.MyPATH+'Results_storage/OIII/'+ID+'_'+sp_binning+'_spaxel_fit.fits', overwrite=True)
     
        
    return storage

import Fitting_tools as emfit

 

def flux_measure_ind(out,wave, mode, use='tot', error=0):
        
    if (mode =='OIII') &( use=='tot'):
        
        fl =( out.eval_components(x= wave)['o3rw_'] + out.eval_components(x= wave)['o3rn_']) *1.333
        
    if (mode =='OIII')&( use=='nar'):
        
        fl = out.eval_components(x= wave)['o3rn_'] *1.333
        
    if (mode =='OIII')&( use=='bro'):
                
        fl = out.eval_components(x= wave)['o3rw_'] *1.333
        
           
    elif (mode =='OIIIs'):
        
        fl = out.eval_components(x= wave)['o3r_']*1.333
        
    
    elif (mode =='H') & (use=='tot'):       
        wvs = np.where( (wave< out.params['Haw_center'].value + 1.5*out.params['Haw_fwhm'].value) & (wave> out.params['Haw_center'].value - 1.5*out.params['Haw_fwhm'].value))[0]
        wave = wave[wvs]
        
        fl = out.eval_components(x= wave)['Han_'] + out.eval_components(x= wave)['Haw_']
    
    
    elif (mode =='H') & (use=='BPT'):
        wvs = np.where( (wave< out.params['Han_center'].value + 1.5*out.params['Han_fwhm'].value) & (wave> out.params['Han_center'].value - 1.5*out.params['Han_fwhm'].value))[0]
        wave = wave[wvs]
        
        fl = out.eval_components(x= wave)['Han_']
        
    elif (mode =='Hb') & (use=='tot'):
        
        fl = out.eval_components(x= wave)['Hbn_'] + out.eval_components(x= wave)['Hbw_']
    
    elif (mode =='Hb') & (use=='BPT'):    
        
        wvs = np.where( (wave< out.params['Hbn_center'].value + 1.5*out.params['Hbn_fwhm'].value) & (wave> out.params['Hbn_center'].value - 1.5*out.params['Hbn_fwhm'].value))[0]
        wave = wave[wvs]
        
        fl = out.eval_components(x= wave)['Hbn_']
        
    
    elif (mode =='Hb') & (use=='nar'):    
        
            
        wvs = np.where( (wave< out.params['Hbn_center'].value + 1.5*out.params['Hbn_fwhm'].value) & (wave> out.params['Hbn_center'].value - 1.5*out.params['Hbn_fwhm'].value))[0]
        wave = wave[wvs]
            
        fl = out.eval_components(x= wave)['Hbn_']
        
    elif (mode =='Hb') & (use=='bro'):    
        
            
        #wvs = np.where( (wave< out.params['Hbw_center'].value + 2*out.params['Hbw_fwhm'].value) & (wave> out.params['Hbw_center'].value - 2*out.params['Hbw_fwhm'].value))[0]
        #wave = wave[wvs]
            
        fl = out.eval_components(x= wave)['Hbw_']
        
    elif (mode =='Hbs'):       
        wvs = np.where( (wave< out.params['Hb_center'].value + 1.5*out.params['Hb_fwhm'].value) & (wave> out.params['Hb_center'].value - 1.5*out.params['Hb_fwhm'].value))[0]
        wave = wave[wvs]
        #    
        fl = out.eval_components(x= wave)['Hb_']    
        
        
    elif (mode =='NII'):       
        wvs = np.where( (wave< out.params['Nr_center'].value + 1.5*out.params['Nr_fwhm'].value) & (wave> out.params['Nr_center'].value - 1.5*out.params['Nr_fwhm'].value))[0]
        wave = wave[wvs] 
        
        fl = out.eval_components(x= wave)['Nr_']*1.333
        
        
    elif (mode =='H') & (use=='broad'):
        #wvs = np.where( (wave< out.params['Haw_center'].value + 1.5*out.params['Haw_fwhm'].value) & (wave> out.params['Haw_center'].value - 1.5*out.params['Haw_fwhm'].value))[0]
        #wave = wave[wvs]
        
        fl = out.eval_components(x= wave)['Haw_']

    import scipy.integrate as scpi
    
    Fit = scpi.simps(fl, wave) * 1e-13
    
    er = np.zeros_like(wave)
    
    if error == 0:
        return Fit
    else:   
        try:
            
            for i in range(len(wave)):
                er[i] =  out.eval_uncertainty(x=wave[i],sigma=1)[0]
            
            Error = Fit - scpi.simps(fl-er, wave)* 1e-13
        
        except:
            Error = 0.1*Fit
            
       
        return Fit, Error

def flux_measure(storage, mode, use='tot'):
    
    wave = storage['obs_wave']
    
    
    if (mode =='OIII'):
        out = storage['1D_fit_OIII_mul']
        
        fl = out.eval_components(x= wave)['o3rw_'] + out.eval_components(x= wave)['o3rn_']
        #fl = fl*1.333
        
    
    elif (mode =='H') & (use=='tot'):       
        out = storage['1D_fit_Halpha_mul']
        
        fl = out.eval_components(x= wave)['Han_'] + out.eval_components(x= wave)['Haw_']
    
    elif (mode =='H') & (use=='BPT'):       
        out = storage['1D_fit_Halpha_mul']
        
        fl = out.eval_components(x= wave)['Han_']
        
    elif (mode =='Hb') & (use=='tot'):
        
        out = storage['1D_fit_Hbeta_mul']
        
        fl = out.eval_components(x= wave)['Hbn_'] + out.eval_components(x= wave)['Hbw_']
    
    elif (mode =='Hb') & (use=='BPT'):        
        out = storage['1D_fit_Hbeta_mul']
        
        fl = out.eval_components(x= wave)['Hbn_']
        
    elif (mode =='NII'):       
        out = storage['1D_fit_Halpha_mul']
        
        fl = out.eval_components(x= wave)['Nr_']

    import scipy.integrate as scpi
    Fit = scpi.simps(fl, wave) * 1e-13
    
    return Fit
    

    

def STD_calc(wave,flux, mode):
    ''' 
    Calculating the STD of a spectrum 
    Parameters
    ----------
    wave : TYPE
        DESCRIPTION.
    flux : TYPE
        DESCRIPTION.
    mode : TYPE
        DESCRIPTION.

    Returns
    -------
    error : TYPE
        DESCRIPTION.
        

    '''
    
    fl = flux.data
    msk = flux.mask    
    em_msk = msk.copy()
    
    if mode=='H':
        cent = 6565.
        side = 400.
        em_line = np.where((wave>6500)&(wave<6650))[0]
        #em_line = np.append(em_line, np.where((wave>6710)&(wave<6720)[0]))
        
        em_msk[em_line] = True
        
        
    elif mode=='OIII':
        cent = 5008.
        side = 150.
        em_line = np.where((wave>4950)&(wave<5030))[0]
        em_line = np.append(em_line, np.where((wave>4840)&(wave<4880))[0])
        
        em_msk[em_line] = True
        
        
        
    elif mode=='Hb':
        cent = 4860.
        side = 100.
        em_line = np.where((wave>4950)&(wave<5030))[0]
        em_line = np.append(em_line, np.where((wave>4840)&(wave<4880))[0])
        
        em_msk[em_line] = True
    
    elif mode=='Hbs':
        cent = 4860.
        side = 100.
        em_line = np.where((wave>4950)&(wave<5030))[0]
        em_line = np.append(em_line, np.where((wave>4840)&(wave<4880))[0])
        
        em_msk[em_line] = True
    
    flx = np.ma.array(data = fl, mask=em_msk)    
    flx = flx[np.where((wave> (cent-side))&(wave< (cent+side)))[0]]
    wv = wave[np.where((wave> (cent-side))&(wave< (cent+side)))[0]]
    
    error = np.ma.std(flx)
    
   
    return error
    
 
def gauss(x,k,mu,sig):
    expo= -((x-mu)**2)/(sig*sig)
    
    y= k* e**expo
    
    return y
       
def Upper_lim_calc(wave,flux,error,z, mode, width_test, fit_u):   
        
    if (mode =='Hb') :        
        ffit = fit_u.eval_components(x= wave)['Hbn_'] + fit_u.eval_components(x= wave)['Hbw_']
        cent = 4865.*(1+z)/1e4
        width_test = width_test/3e5* (4862.*(1+z)/1e4)/2.35 #2.35 is to convert between FWHM and sigma
        
        
    
    if (mode =='Hbs') :        
        ffit = fit_u.eval_components(x= wave)['Hb_'] + fit_u.eval_components(x= wave)['Hb_']
        cent = 4865.*(1+z)/1e4
        width_test = width_test/3e5* (4862.*(1+z)/1e4)/2.35 #2.35 is to convert between FWHM and sigma
        
        
        
    
    flux_sub  = flux- ffit
 
    N_at = 40
    peak_test = np.linspace(1.5, 10, N_at)* error.data[0]    

    for i in range(N_at):
        flux_test = flux_sub + gauss(wave, peak_test[i], cent, width_test)

        if mode=='Hb':
            New,c = emfit.fitting_Hbeta_mul(wave, flux_test, error,z, decompose= np.array([cent, 10, 0, width_test]))                                  
            SNR,chi2 = SNR_calc(flux_test, wave, New, 'Hb',z, mul=1)
            #print 'SNR is ' ,SNR
        
        if mode=='Hbs':
            New,c = emfit.fitting_Hbeta_mul(wave, flux_test, error,z, decompose= np.array([cent, 10, 0, width_test]))                                  
            SNR,chi2 = SNR_calc(flux_test, wave, New, 'Hb',z, mul=1)
            #print 'SNR is ' ,SNR
        
        if SNR>3.:
            break
    f,(ax2) = plt.subplots(1)
 
    ax2.plot(wave/(1+z)*1e4, flux_test.data, color='grey',drawstyle='steps-mid', label='Unmasked')
    
    #Substracted
    ax2.plot(wave[np.invert(flux_sub.mask)]/(1+z)*1e4, flux_sub[np.invert(flux_sub.mask)], color='blue',drawstyle='steps-mid', label='Substracted')
    
    # Added
    ax2.plot(wave[np.invert(flux_test.mask)]/(1+z)*1e4, flux_test[np.invert(flux_test.mask)], color='orange',drawstyle='steps-mid', label='Upper lim')
    
    #Original
    ax2.plot(wave[np.invert(flux.mask)]/(1+z)*1e4, flux[np.invert(flux.mask)], color='red',drawstyle='steps-mid', label='Original')
    
    ax2.set_xlim(4862-40, 4862+40)
    ax2.set_ylim(-0.1, 0.5)
    ax2.legend(loc='best')
    
    if mode =='Hb':
        t=1
        #emplot.plotting_Hbeta(wave, flux_test, ax2, New,z, 'mul', title=1)
    if mode=='Hbs':
        Flux_up = flux_measure_ind(New,wave, 'Hb' , use='BPT')
    else:
        Flux_up = flux_measure_ind(New,wave, mode , use='BPT')
        
    print ('Flux upper limit is ',Flux_up)
    
    return Flux_up
 


def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx   

def W80_mes(out, mode, plot):
    import scipy.integrate as scpi
    
    if mode=='OIII':
        N= 5000
        
        cent = out.params['o3rn_center'].value
        
        bound1 =  cent + 4000/3e5*cent
        bound2 =  cent - 4000/3e5*cent
              
        wvs = np.linspace(bound2, bound1, N)
        
        try:            
            y = out.eval_components(x=wvs)['o3rw_'] + out.eval_components(x=wvs)['o3rn_']
        
        except:
            y = out.eval_components(x=wvs)['o3rn_']
                            
        Int = np.zeros(N-1)
        
        for i in range(N-1):
            
            Int[i] = scpi.simps(y[:i+1], wvs[:i+1]) * 1e-13
                    
        Int = Int/max(Int)

        ind10 = np.where( (Int>0.1*0.992)& (Int<0.1/0.992) )[0][0]            
        ind90 = np.where( (Int>0.9*0.995)& (Int<0.9/0.995) )[0][0] 
        ind50 = np.where( (Int>0.5*0.992)& (Int<0.5/0.992) )[0][0]            
        
        wv10 = wvs[ind10]
        wv90 = wvs[ind90]
        wv50 = wvs[ind50]
        
        v10 = (wv10-cent)/cent*3e5
        v90 = (wv90-cent)/cent*3e5
        v50 = (wv50-cent)/cent*3e5
        
        
        w80 = v90-v10
        
        
        if plot==1:
            f,ax1 = plt.subplots(1)
            ax1.plot(wvs,y, 'k--')
        
        
            g, ax2 = plt.subplots(1)
            ax2.plot(wvs[:-1], Int)
        
            ax2.plot(np.array([bound2,bound1]), np.array([0.9,0.9]), 'r--')
            ax2.plot(np.array([bound2,bound1]), np.array([0.1,0.1]), 'r--')
        
            ax2.plot(np.array([cent,cent]), np.array([0,1]), 'b--')
        
            ax1.plot(np.array([wv10,wv10]), np.array([0, max(y)]), 'r--')
            ax1.plot(np.array([wv90,wv90]), np.array([0, max(y)]), 'r--')
            
    
    elif mode=='CO':
        N= 5000
        
        cent = out.params['COn_center'].value
        
        bound1 =  cent + 4000/3e5*cent
        bound2 =  cent - 4000/3e5*cent
              
        wvs = np.linspace(bound2, bound1, N)
        
        try:            
            y = out.eval_components(x=wvs)['COn_'] + out.eval_components(x=wvs)['COb_']
        
        except:
            y = out.eval_components(x=wvs)['COn_']
                            
        Int = np.zeros(N-1)
        
        for i in range(N-1):
            
            Int[i] = scpi.simps(y[:i+1], wvs[:i+1]) * 1e-13
                    
        Int = Int/max(Int)
        try:
            
            ind10 = find_nearest(Int, 0.1)          
               
        
        except:
            print( np.where( (Int>0.1*0.991)& (Int<0.1/0.991) )[0]  )
            plt.figure()
            plt.plot(wvs[:-1], Int)
           
            
            
        ind90 = np.where( (Int>0.9*0.995)& (Int<0.9/0.995) )[0][0] 
        ind50 = np.where( (Int>0.5*0.992)& (Int<0.5/0.992) )[0][0] 
        
        wv10 = wvs[ind10]
        wv90 = wvs[ind90]
        wv50 = wvs[ind50]
        
        v10 = (wv10-cent)/cent*3e5
        v90 = (wv90-cent)/cent*3e5
        v50 = (wv50-cent)/cent*3e5
        
        
        w80 = v90-v10
        
        
        if plot==1:
            f,ax1 = plt.subplots(1)
            ax1.plot(wvs,y, 'k--')
        
        
            g, ax2 = plt.subplots(1)
            ax2.plot(wvs[:-1], Int)
        
            ax2.plot(np.array([bound2,bound1]), np.array([0.9,0.9]), 'r--')
            ax2.plot(np.array([bound2,bound1]), np.array([0.1,0.1]), 'r--')
        
            ax2.plot(np.array([cent,cent]), np.array([0,1]), 'b--')
        
            ax1.plot(np.array([wv10,wv10]), np.array([0, max(y)]), 'r--')
            ax1.plot(np.array([wv90,wv90]), np.array([0, max(y)]), 'r--')
    
    return v10,v90, w80, v50


def W68_mes(out, mode, plot):
    import scipy.integrate as scpi
    
    if mode=='OIII':
        N= 5000
        
        cent = out.params['o3rn_center'].value
        
        bound1 =  cent + 4000/3e5*cent
        bound2 =  cent - 4000/3e5*cent
              
        wvs = np.linspace(bound2, bound1, N)
        
        try:            
            y = out.eval_components(x=wvs)['o3rw_'] + out.eval_components(x=wvs)['o3rn_']
        
        except:
            y = out.eval_components(x=wvs)['o3rn_']
                            
        Int = np.zeros(N-1)
        
        for i in range(N-1):
            
            Int[i] = scpi.simps(y[:i+1], wvs[:i+1]) * 1e-13
                    
        Int = Int/max(Int)

        ind16 = np.where( (Int>0.16*0.992)& (Int<0.16/0.992) )[0][0]            
        ind84 = np.where( (Int>0.84*0.995)& (Int<0.84/0.995) )[0][0] 
        ind50 = np.where( (Int>0.5*0.992)& (Int<0.5/0.992) )[0][0]            
        
        wv10 = wvs[ind16]
        wv90 = wvs[ind84]
        wv50 = wvs[ind50]
        
        v10 = (wv10-cent)/cent*3e5
        v90 = (wv90-cent)/cent*3e5
        v50 = (wv50-cent)/cent*3e5
        
        
        w80 = v90-v10
        
        
        if plot==1:
            f,ax1 = plt.subplots(1)
            ax1.plot(wvs,y, 'k--')
        
        
            g, ax2 = plt.subplots(1)
            ax2.plot(wvs[:-1], Int)
        
            ax2.plot(np.array([bound2,bound1]), np.array([0.9,0.9]), 'r--')
            ax2.plot(np.array([bound2,bound1]), np.array([0.1,0.1]), 'r--')
        
            ax2.plot(np.array([cent,cent]), np.array([0,1]), 'b--')
        
            ax1.plot(np.array([wv10,wv10]), np.array([0, max(y)]), 'r--')
            ax1.plot(np.array([wv90,wv90]), np.array([0, max(y)]), 'r--')
            
    
    elif mode=='CO':
        N= 5000
        
        cent = out.params['COn_center'].value
        
        bound1 =  cent + 4000/3e5*cent
        bound2 =  cent - 4000/3e5*cent
              
        wvs = np.linspace(bound2, bound1, N)
        
        try:            
            y = out.eval_components(x=wvs)['COn_'] + out.eval_components(x=wvs)['COb_']
        
        except:
            y = out.eval_components(x=wvs)['COn_']
                            
        Int = np.zeros(N-1)
        
        for i in range(N-1):
            
            Int[i] = scpi.simps(y[:i+1], wvs[:i+1]) * 1e-13
                    
        Int = Int/max(Int)
        try:
            
            ind10 = find_nearest(Int, 0.16)          
               
        
        except:
            print( np.where( (Int>0.16*0.991)& (Int<0.16/0.991) )[0]  )
            plt.figure()
            plt.plot(wvs[:-1], Int)
           
            
            
        ind90 = np.where( (Int>0.84*0.995)& (Int<0.84/0.995) )[0][0] 
        ind50 = np.where( (Int>0.5*0.992)& (Int<0.5/0.992) )[0][0] 
        
        wv10 = wvs[ind10]
        wv90 = wvs[ind90]
        wv50 = wvs[ind50]
        
        v10 = (wv10-cent)/cent*3e5
        v90 = (wv90-cent)/cent*3e5
        v50 = (wv50-cent)/cent*3e5
        
        
        w80 = v90-v10
        
        
        if plot==1:
            f,ax1 = plt.subplots(1)
            ax1.plot(wvs,y, 'k--')
        
        
            g, ax2 = plt.subplots(1)
            ax2.plot(wvs[:-1], Int)
        
            ax2.plot(np.array([bound2,bound1]), np.array([0.9,0.9]), 'r--')
            ax2.plot(np.array([bound2,bound1]), np.array([0.1,0.1]), 'r--')
        
            ax2.plot(np.array([cent,cent]), np.array([0,1]), 'b--')
        
            ax1.plot(np.array([wv10,wv10]), np.array([0, max(y)]), 'r--')
            ax1.plot(np.array([wv90,wv90]), np.array([0, max(y)]), 'r--')
    
    return v10,v90, w80, v50


def Av_calc(Falpha, Fbeta):
    '''Calculating Av from Halpha flux and Hbeta flux
    '''
    
    Av = 1.964*4.12*np.log10(Falpha/Fbeta/2.86)
    
    return Av

def Flux_cor(Flux, Av, lam= 0.6563):
    '''Correcting Halpha flux based on Av
    '''
    
    Ebv = Av/4.12
    Ahal = 3.325*Ebv
    
    F = Flux*10**(0.4*Ahal)
    
    return F

def Sub_QSO(storage_H):
    ID = storage_H['X-ray ID']
    z = storage_H['z_guess']
    flux = storage_H['flux'].copy()
    Mask = storage_H['sky_clipped']
    Mask_1D = storage_H['sky_clipped_1D']
    out = storage_H['1D_fit_Halpha_mul']
    
    wv_obs = storage_H['obs_wave'].copy()
    
    rst_w = storage_H['obs_wave'].copy()*1e4/(1+storage_H['z_guess'])
        
    center =  storage_H['Median_stack_white_Center_data'][1:3].copy()
    sig = storage_H['Signal_mask'][0,:,:].copy()
    sig[:,:] = True
    
    shapes = storage_H['dim']
    # This choose spaxel within certain radius. Then sets it to False since we dont mask those pixels
    for ix in range(shapes[0]):
        for iy in range(shapes[1]):
            dist = np.sqrt((ix- center[1])**2+ (iy- center[0])**2)
            if dist< 6:
                sig[ix,iy] = False
        
    Mask_em = Mask_1D.copy()
    Mask_em[:] = False
        
    wmin = out.params['Nb_center'].value - 1*out.params['Nb_sigma'].value
    wmax = out.params['Nr_center'].value + 1*out.params['Nr_sigma'].value

    
        
    em_nar = np.where((wv_obs>wmin) & (wv_obs<wmax))[0]        
    Mask_em[em_nar] = True
    
    
    
    wmin = (6718.-5)*(1+z)/1e4
    wmax = (6732.+5)*(1+z)/1e4
    
    em_nar = np.where((wv_obs>wmin) & (wv_obs<wmax))[0]        
    Mask_em[em_nar] = True

       
    comb = np.logical_or(Mask_1D, Mask_em)
    
    shapes = storage_H['dim']
    Cube = np.zeros((shapes[2], shapes[0], shapes[1]))
    C_ms = Cube[:,:,:].copy()
    C_ms[:,:,:] = False
    Cube = np.ma.array(data=Cube, mask = C_ms)
    
    BLR_map = np.zeros((shapes[0], shapes[1]))
    
    
    
    x = range(shapes[0])
    y = range(shapes[1])
    
    #x = np.linspace(30,60,31)
    #y = np.linspace(30,60,31)
    
    x = np.array(x, dtype=int)
    y = np.array(y, dtype=int)
        
    ls,ax = plt.subplots(1)
    pdf_plots = PdfPages(ph.MyPATH+'Graphs/'+ID+'_SUB_QSO.pdf') 
    #import progressbar
    for i in x:# progressbar.progressbar(x):
        print (i,'/',x[-1])
        for j in y:
            flx_pl = flux[:,i,j]
            flxm = np.ma.array(data=flx_pl, mask=comb)
            
                
            error = STD_calc(wv_obs/(1+z)*1e4,flxm, 'H')* np.ones(len(flxm))
            
            plt.figure
            
            plt.plot(wv_obs, flxm.data)
            
            try:
                
            
                out = emfit.sub_QSO(wv_obs, flxm, error,z, storage_H['1D_fit_Halpha_mul'])            
                sums = flxm.data-(out.eval_components(x=wv_obs)['Ha_']+out.eval_components(x=wv_obs)['linear'])
                Cube[:,i,j] = sums
                
                BLR_map[i,j] = np.sum(out.eval_components(x=wv_obs)['Ha_'])
                    
                    
                if sig[i,j]==False:
                    ax.plot(wv_obs, flxm.data, label='data')
                    ax.plot(wv_obs,out.eval_components(x=wv_obs)['Ha_']+out.eval_components(x=wv_obs)['linear'], label='model')
                    ax.plot(wv_obs, sums, label='Res')
                        
                    ax.legend(loc='best')
                    wmin = out.params['Ha_center'].value - 3*out.params['Ha_sigma'].value
                    wmax = out.params['Ha_center'].value + 3*out.params['Ha_sigma'].value
                        
                    ax.set_xlim(wmin,wmax)
                        
                    pdf_plots.savefig()
                    
                    ax.clear()
            except:
                Cube[:,i,j] = flxm.data
                
                BLR_map[i,j] =0
                
                print (i,j, ' BLR sub fail')
                
                
                
    pdf_plots.close()

    Cube_ready = np.ma.array(data= Cube, mask= storage_H['flux'].copy().mask)               
    storage_new = storage_H.copy()   
    storage_new['flux'] = Cube_ready
    storage_new['BLR_map'] = BLR_map
        
        
    prhdr = storage_H['header']
    hdu = fits.PrimaryHDU(Cube_ready.data, header=prhdr)
    hdulist = fits.HDUList([hdu])
    
    hdulist.writeto(ph.MyPATH+'Sub_QSO/'+ID+'.fits', overwrite=True)
    
    return storage_new




