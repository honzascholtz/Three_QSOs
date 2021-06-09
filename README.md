# The impact of ionised outflows from z~2.5 quasars is not through instantaneous in-situ quenching: the evidence from ALMA and VLT/SINFONI
This is code to analyze the Halpha, [OIII] and FIR emission of three QSOs - 2QZJ002830.4-2817, LBQS0109+0213 and HB89 0329-385 from work:

The impact of ionised outflows from z$\sim$2.5 quasars is not through instantaneous in-situ quenching: the evidence from ALMA and VLT/SINFONI

by:

J. Scholtz,  C.M. Harrison, D.J. Rosario, D.M. Alexander, K.K. Knudsen, F. Stanley, Chian-Chou Chen, D. Kakkad, V. Mainieri and J. Mullaney

DOI: 10.1093/mnras/stab1631

arXiv: 

All the data and products are available to download at https://drive.google.com/drive/folders/1K2Fnd-H_8nu-SKr84noMplyia3TXRl-B?usp=sharing. 

Please download the whole package and change the MyPATH in Tools_path.py to point to the folder containing the code and Works_folder. 

I have included all the by-products so you can run any scripts and it should generate the figure in the paper. If you would like to start from scratch, please run the script in the order they are named. This will allow to create all the products necessary to plot the paper figures. 

This code has been tested on Python 3.7, numpy 1.15.3, scipy 1.2.1. Preliminary testing suggested that the code will NOT work with numpy 1.16.0 or latee.


#Description of the products:

# ALMA data

Work_folder/ALMA/Final_set/new/* - Imaged ALMA band 7 continuum data at natural resolution

Work_folder/ALMA/UV_data/* - Collapsed uv-visibiities for each object after subtracting any other sources in the field


# CATALOGUES 

Work_folder/Catalogues/Chains/* - Result of the FortesFit MCMC chains

Work_folder/Catalogues/3QSO_phot_convert.fits - Photometry catalogue used for SED fitting

Work_folder/Catalogues/QSOs.fit_results.fits - Results of the SED fitting for these QSOs

Work_folder/Catalogues/Schulze19_QSOs.fit_results - SED fitting of targets from Schulze+19

Work_folder/Catalogues/Stanley_data.fits - Results from Stacking from Stanley+17

Work_folder/Catalogues/KASHZ_OIII - OIII results data from KASHz survey (Harrison+16)

Work_folder/Catalogues/Scholtz_2020_OIII - OIII results from Scholtz+2020

Work_folder/Catalogues/QSO_netzer - QSO results from Netzer+04

Work_folder/Catalogues/QSO_shen - QSO results from Shen+04

Work_folder/Catalogues/ALPAKA_XMM - results from the ALPAKA survey cross matched with XMM catalogue

# SINFONI DATA
Work_folder/SINFONI/* H & K SINFONI cubes

Work_folder/Sub_QSO/* Cubes with subtracted BLR and continuum components as described in the Appendix 

# RESULTS_storage - Code products

Work_folder/Results_storage/Growth/* Curves-of-growth results stored

Work_folder/Results_storage/Halpha/* .fits - results of the spaxel-by-spaxel fitting 

Work_folder/Results_storage/Halpha/* _res.fits - Residual cube from the spaxel-by-spaxel fitting 

Work_folder/Results_storage/Halpha/* _mod.fits - Model from the spaxel-by-spaxel fitting 

Work_folder/Results_storage/imfit/* Results from the CASA's imfit fitting

Work_folder/Results_storage/OIII/* Results of the OIII spaxel-by-spaxel fitting

Work_folder/Results_storage/Properties/* The storage of dictionaries used the code

Work_folder/Results_storage/Spectrums/* Storage of the different inner and outer spectra

Work_folder/Results_storage/Spectrums/Regions/* Storage of the regional spectra for [OIII] and Halpha

Work_folder/Results_storage/Sub_qso_map_hal/* Map of the narrow Halpha emission after collapsing the residual cube from subtracting the BLR and continuum only

Work_folder/Results_storage/* PSF_sub_smt.txt - Narrow Halpha map with PSF subtracted

# GRAPHS

Work_folder/Results_storage/Graphs/
Work_folder/Results_storage/Graphs/








