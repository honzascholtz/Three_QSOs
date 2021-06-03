import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import emcee
import pdb
import corner
import math
import sys

from scipy import special as sp
from tqdm import tqdm

import warnings
warnings.filterwarnings("ignore")

# Model for the log-normal
def modelPDF(x, mu, sigma):
	pdf = 1./sigma/np.sqrt(2.*3.14) * np.exp(-((x)-mu)**2./2./sigma**2.)
	return pdf

# Model for the log-normal for upper and lower limits
def modelCDF(x, mu, sigma):
	cdf = 0.5 + 0.5 * sp.erf(((x)-mu)/np.sqrt(2.)/sigma)
	return cdf

# Prior distributions
def lnprior(theta, Pmu, Psig):
    mu,sigma = theta

    if  Pmu[0] < mu < Pmu[1] and Psig[0] < sigma < Psig[1]:
    	return 0.0
    return -np.inf

def lnprob(theta, x, ex, Pmu, Psig):
	lp = lnprior(theta, Pmu, Psig)
	if not np.isfinite(lp):
	    return -np.inf
	return lp + lnlike(theta, x, ex, Pmu, Psig)

# Likelihood function
def lnlike(theta, x, ex, Pmu, Psig):
	mu = theta[0]
	sigma = theta[1]
	
	o = np.where((ex >= 0) & (ex != 99.))[0]
	likdet = np.sum(np.log(modelPDF(x[o],mu,sigma)))
	
	o = np.where(ex < 0)[0]
	likundet = np.sum(np.log(modelCDF(abs(x[o]),mu,sigma)))

	o = np.where(ex == 99.)[0]
	liklowlim = np.sum(np.log(1. - modelCDF(abs(x[o]),mu,sigma)))

	return likdet + likundet + liklowlim


# Name of the file in csv format
# Values need to be non-negative
# Upper limits need to be -99. in the uncertainties (the value in the parameter column being the upper limit)
# Lower limits need to be 99. in the uncertainties(the value in the parameter column being the lower limit)
# FOR TEST df = pd.read_csv('ptest.csv')
df = pd.read_csv('/Users/jansen/Google Drive/Astro/Four_Quasars/Catalogues/Schulze_dist.csv')

# Name of the column that contains the parameter (can be a list)
# The name of the uncerainties need to be the name of the parameter preceeded by an e
# (e.g. if parameter is SFR the columns containing the uncertainties needs to be called ESFR)
# FOR TEST name_param = ['SFR']
name_param = ['SFR']

# Number of Monte Carlo for the uncertainties
Nmc_unc = 100
# Number of Monte Carlo for the sampling
Nmc = 5000

#Values of the priors on the mu and sigma
Pmu = [-4., 4.]
Psig = [0.01, 3.]

for name_i in name_param:

	try:
		del flat_chain
	except:
		pass

	for Nmc_i in tqdm(range(0, Nmc_unc)):

		try:
			# Extract the parameters (wether being detected or not)
			p = df[name_i].values
			# Extract the uncertainties (wether being detected or not)
			ep = df['E'+name_i].values
		except:
			print ("_________________________")
			print ('Make sure that you have no space in between the names in the header of your csv file '\
				+ 'i.e. "SFR,ESFR" if ok while "SFR,   ESFR"   is not ok (hit ctr+c and enter once you read this)')
			print( "_________________________")
			print( "_________________________")
			print( "_________________________")
			print( "_________________________")
			pdb.set_trace()

		# Sample from the observations
		psamp = p.copy()
		# Select only the detected values
		o = np.where((ep != -99.) & (ep != 99.))[0]
		# Randomly draw values within the uncertainties
		psamp[o] = np.random.normal(p[o], ep[o],len(o))
		# If it happens that the new value is negative, replace it with an upper limit
		o = np.where(psamp<0.)[0]
		if len(o) > 0.:
			psamp[o] = p[o]
			ep[o] = -99.

		# MCMC
		ndim, nwalkers = 2, 20
		p0 = np.zeros(shape=(nwalkers, ndim))
		p0[:,0] = np.random.uniform(low = Pmu[0], high = Pmu[1], size=nwalkers)
		p0[:,1] = np.random.uniform(low = Psig[0], high = Psig[1], size=nwalkers)

		sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(psamp, ep, Pmu, Psig))

		# Run the MCMC
		chain = sampler.run_mcmc(p0, Nmc)

		#Result of the fit
		try:
			flat_chain = np.concatenate((flat_chain, sampler.chain[:,int(0.1*Nmc):,:].reshape([-1,ndim])), axis = 0)
		except:
			flat_chain = sampler.chain[:,int(0.1*Nmc):,:].reshape([-1,ndim])

	mu = np.median(flat_chain[:,0])
	mu_err = np.std(flat_chain[:,0])
	sigma = np.median(flat_chain[:,1])
	sigma_err = np.std(flat_chain[:,1])

	#Corner plot of the fit
	fig = corner.corner(flat_chain, labels=[r"$\mu$", r"$\sigma$"], truths = [mu,sigma])

	# Save the Corner plot (just put the PATH)
	# FOR TEST fig.savefig('./test.pdf', dpi=300)
	fig.savefig('/Users/jansen/Google Drive/Astro/Four_Quasars/Graphs/Schulze_CORNER.pdf', dpi=300)

	# Save the chain
	df_chain = pd.DataFrame(flat_chain)
	df_chain.columns = ['mu','sigma']
	
	# Save the chain as a csv file (just change the PATH)
	#FOR TEST df_chain.to_csv('./test.csv', index = 0, columns = ['mu','sigma'])
	df_chain.to_csv('/Users/jansen/Google Drive/Astro/Four_Quasars/Schulze_CHAIN.csv', index = 0, columns = ['mu','sigma'])

from astropy.table import Table, join, vstack
Control_res = Table.read('/Users/jansen/Google Drive/Astro/Four_Quasars/Schulze_CHAIN.csv', format='csv')



Cmu_ar = Control_res['mu']
Csig_ar = Control_res['sigma']

Cmu = np.median(Cmu_ar)
Csig = np.median(Csig_ar)