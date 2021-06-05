#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 23 10:44:18 2021

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

# =============================================================================
# This was read from the PDF files and saved in those csv files as suggested by David
# =============================================================================
q_2QZJ = Table.read(ph.MyPATH+'2qzj.csv', format='csv')

Rat = q_2QZJ[3][1]/q_2QZJ[0][1]*100
Rat_up = q_2QZJ[2][1]/q_2QZJ[1][1]*100
print(Rat, Rat_up)

q_LBQS = Table.read(ph.MyPATH+'LBQS.csv', format='csv')

Rat = q_LBQS[3][1]/q_LBQS[0][1]*100
Rat_up = q_LBQS[2][1]/q_LBQS[1][1]*100
print(Rat, Rat_up)


# =============================================================================
# Creating class to read FortesFit results chains - taken directly from FortesFit source code
# =============================================================================
import h5py

class FortesFitResult:
	""" Representation of the FORTES-FIT MCMC output for a single object """
	
	
	def __init__(self,FortesFit_OutFile, BurnIn = 10000, old=False):
		""" Read in the FORTES-FIT MCMC outputs for a single object 
			
			FortesFit_OutFile: The output HDF5 file of FortesFit_FitSingle_emcee
			BurnIn: number of initial samples to exclude to allow for convergence of the chains. Default = 10000 
			
		"""
				
		FitFile = h5py.File(FortesFit_OutFile, 'r')	 #  Open the HDF5 output file	

		self.objectname = FitFile.attrs['Object Name'].decode()
		self.fit_description = FitFile.attrs['Description'].decode()
			
		priors = OrderedDict()  #  A dictionary that stores priors for all parameters, including redshift
			
		# Object specific attributes		
		priors.update({'Redshift':FitFile.attrs['Redshift']})   

		self.fit_filterids = FitFile['Photometry/FilterIDs'][()] # Filters of the SED used in the fit
		self.fit_fluxes = FitFile['Photometry/Fluxes'][()] # Filters of the SED used in the fit
		self.fit_fluxerrors = FitFile['Photometry/FluxErrors'][()] # Errors on the fluxes in erg/s/cm^2/micron

		# Basic attributes about the overall fits, from the main file metadata
		self.fit_modelids  = FitFile['Model'].attrs['ModelIDs']  # Models used to fit the SED

		for modelid in self.fit_modelids:
			subgroupname = 'Model/Model{0:2d}'.format(modelid)
			for uparam in FitFile[subgroupname].keys():
				dataname = subgroupname+'/'+uparam
				priors.update({uparam:FitFile[dataname][()]})
		
		self.priors = priors
		paramnames = np.array(list(priors.keys()))	
		
		self.fit_parameter_names  = np.core.defchararray.decode(FitFile['Chain/Varying_parameters'][()]) # Ordered list of parameters that were fit
		matchindex = np.zeros(len(self.fit_parameter_names),dtype='i2')
		for iparam,param in enumerate(self.fit_parameter_names):
			index, = np.where(paramnames == param)
			matchindex[iparam] = index[0]
		self.fit_parameter_indices = matchindex
		
		self.burn_in = BurnIn

		if old:
			tempchains  = FitFile['Chain/emcee_chain'][()]
			self.chains = tempchains.reshape((tempchains.shape[0]*tempchains.shape[1],tempchains.shape[2]),order='F')
#			self.chains = tempchains.reshape((-1,len(self.fit_parameter_names)))  #  Store the entire EMCEE output chain
			self.all_samples = self.chains[BurnIn:,:]
		else:	
			self.chains = FitFile['Chain/posterior_chain'][()]  #  Store the entire posterior output chain
			self.all_samples = self.chains[BurnIn:,:]
		# Future: add a warning here in case chains are too short for reasonable statistics of posterior PDF
				
		FitFile.close()  #  Close the HDF5 output file
		
		# Use marginalised posterior PDF to calculate Kullback-Leibler Divergence wrt priors of fitted parameters
		posteriors = OrderedDict()  #  A dictionary that stores marginalised posteriors for all parameters, including redshift
		self.fit_parameter_KLD  = np.zeros(len(self.fit_parameter_names))
		for iparam,param in enumerate(self.fit_parameter_names):
			# Use a 30 bin histogram to obtain a marginalised distribution from the posterior samples
			pdfhist = np.histogram(self.all_samples[:,iparam],bins=30)
#			kde = gaussian_kde(pdfhist[0]) # Use KDE to get a smoothed version of the PDF that overcomes zero counts
			post_x = 0.5*(pdfhist[1][0:-1]+pdfhist[1][1:])
#			post_y = kde(post_x)
			post_y = pdfhist[0]
			prior = self.priors[param]
			pk = post_y
			qk = np.interp(post_x,prior[0,:],prior[1,:],left=0.0,right=0.0)
			self.fit_parameter_KLD[iparam] = entropy(pk,qk=qk)
			
#			posteriors.update({param:np.stack([post_x,post_y])})
#			posteriors.update({param:process_PDF(np.stack([post_x,post_y]))})
#		self.marginalised_posteriors = posteriors
			

		# Best-fit or fixed values for each parameter
		perc_pars = {}
		# First variable parameters
		for iparam,param in enumerate(self.fit_parameter_names):
			perc_pars.update({param:(np.percentile(self.all_samples[:,iparam],[50])[0],'Fit')})
		# Then fixed parameters
		for param in self.priors.keys():
			if self.priors[param].shape[0] == 1:
				# Fixed parameters have single element prior arrays
				perc_pars.update({param:(self.priors[param][0],'Fixed')})
		
		self.bestfit_parameters = perc_pars		

		# Redshift of the object
		self.redshift = self.bestfit_parameters['Redshift'][0]
			
	
	def	percentiles(self, Quantiles = [16,50,84]):
		"""	Calculate percentile ranges on varying parameters for a set of input quantiles
	
			Quantiles: list or array of quantiles, Default is equivalent to -1 sigma, median, +1 sigma
		"""
		
		perc_pars = {}
		# Variable parameters only
		for iparam,param in enumerate(self.fit_parameter_names):
			perc_pars.update({param:np.percentile(self.all_samples[:,iparam],Quantiles)})
		
		return perc_pars

# =============================================================================
# Plotting LBQS degeneracies
# =============================================================================
file = ph.MyPATH+'Catalogues/Chains/LBQS.FortesFit_output.multinest.hdf5'



FitFile = h5py.File(file, 'r')

ID = FitFile.attrs['Object Name'].decode()
Info = FitFile.attrs['Description'].decode()
Filters = FitFile['Photometry/FilterIDs'][()]

Fluxes = FitFile['Photometry/Fluxes'][()]
Fluxes_er  = FitFile['Photometry/FluxErrors'][()]

Para_names = np.core.defchararray.decode(FitFile['Chain/Varying_parameters'][()])
print(Para_names)
burnin = 500
Chains = FitFile['Chain/posterior_chain'][()][burnin:,:]



binss = list()
binss.append(np.linspace(44,46.8,10))
binss.append(np.linspace(46.5,47.5,10))


hist_big, x_big,y_big  = np.histogram2d(Chains[:,7],Chains[:,2], bins=binss)
hist_big = hist_big/np.max(hist_big)

f, (ax1,ax2) = plt.subplots(2, sharex=True, sharey=True)
levels_f_plt = (0.68,0.95,0.99)

lens=25
alpha_c = 0.75
width_c = 2.


ax1.contour(y_big[:-1],x_big[:-1], 1-hist_big, 
            levels= levels_f_plt, 
            colors= 'blue', 
            linewidths= width_c, 
            alpha=alpha_c)#,
            #cmap=plt.cm.bone)

ax1.hlines( np.percentile(Chains[:,7],50), 46,48, 'k', linestyle='dashed')
ax1.hlines( np.percentile(Chains[:,7],16), 46,48, 'k', linestyle='dashed',alpha=0.6)
ax1.hlines( np.percentile(Chains[:,7],84), 46,48, 'k', linestyle='dashed',alpha=0.6)


ax1.text(46.52, np.percentile(Chains[:,7],84)+0.05 ,'84th       %.1f' %(np.percentile(Chains[:,7],84)))
ax1.text(46.52, np.percentile(Chains[:,7],50)+0.05 ,'50th       %.1f' %(np.percentile(Chains[:,7],50)))
ax1.text(46.52, np.percentile(Chains[:,7],16)+0.05 ,'16th       %.1f' %(np.percentile(Chains[:,7],16)))


ax1.text(46.59, np.percentile(Chains[:,7],84)+0.05 ,'%')
ax1.text(46.59, np.percentile(Chains[:,7],50)+0.05 ,'%')
ax1.text(46.59, np.percentile(Chains[:,7],16)+0.05 ,'%')



# =============================================================================
# Plotting 2QZJ degenracies. 
# =============================================================================
file = ph.MyPATH+'Catalogues/Chains/2QZJ.FortesFit_output.multinest.hdf5'

FitFile = h5py.File(file, 'r')

ID = FitFile.attrs['Object Name'].decode()
Info = FitFile.attrs['Description'].decode()
Filters = FitFile['Photometry/FilterIDs'][()]

Fluxes = FitFile['Photometry/Fluxes'][()]
Fluxes_er  = FitFile['Photometry/FluxErrors'][()]

Para_names = np.core.defchararray.decode(FitFile['Chain/Varying_parameters'][()])
print(Para_names)
burnin = 500
Chains = FitFile['Chain/posterior_chain'][()][burnin:,:]

hist_big, x_big,y_big  = np.histogram2d(Chains[:,7],Chains[:,2], bins=binss)
hist_big = hist_big/np.max(hist_big)


levels_f_plt = (0.68,0.95,0.99)

lens=25
alpha_c = 0.75
width_c = 2.


ax2.contour(y_big[:-1],x_big[:-1], 1-hist_big, 
            levels= levels_f_plt, 
            colors= 'red', 
            linewidths= width_c, 
            alpha=alpha_c)#,
            #cmap=plt.cm.bone)

ax2.hlines( np.percentile(Chains[:,7],50), 46,48, 'k', linestyle='dashed')
ax2.hlines( np.percentile(Chains[:,7],16), 46,48, 'k', linestyle='dashed', alpha=0.6)
ax2.hlines( np.percentile(Chains[:,7],84), 46,48, 'k', linestyle='dashed', alpha=0.6)

ax2.text(47.2, np.percentile(Chains[:,7],84)+0.05 ,'84th       %.1f' %(np.percentile(Chains[:,7],84)))
ax2.text(47.2, np.percentile(Chains[:,7],50)+0.05 ,'50th       %.1f' %(np.percentile(Chains[:,7],50)))
ax2.text(47.2, np.percentile(Chains[:,7],16)+0.05 ,'16th       %.1f' %(np.percentile(Chains[:,7],16)))


ax2.text(47.27, np.percentile(Chains[:,7],84)+0.05 ,'%')
ax2.text(47.27, np.percentile(Chains[:,7],50)+0.05 ,'%')
ax2.text(47.27, np.percentile(Chains[:,7],16)+0.05 ,'%')

ax1.set_ylabel(r'log$_{\rm 10}$ L$_{\rm IR, SF}$ (erg/s)')
ax2.set_ylabel(r'log$_{\rm 10}$ L$_{\rm IR, SF}$ (erg/s)')
ax2.set_xlabel(r'log$_{\rm 10}$ L$_{\rm IR, AGN}$ (erg/s)')

ax1.tick_params(direction='in')
ax2.tick_params(direction='in')

ax1.set_ylim(44,47)
ax1.set_xlim(46.5,47.45)

plt.savefig(ph.MyPATH+'Graphs/Paper_plots/PDFs_test.pdf', bbox_inches='tight')
    


plt.show()

