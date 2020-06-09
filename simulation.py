import numpy as np
from numpy import random
from prior import prior_dist
from source_functions import *
from plot import plot_samples
import time

def Simulation(ps, par_dict, num_param, redshifts=None, luminosities=None,
obs_pf=None, detector=None, options=None, vol_arr=None, kc=None, dl=None,
plot_GRB=None):

    # Plotting options
    if plot_GRB is None:
        plot_GRB = options.getboolean('plot_GRB')
    plot_func = options.getboolean('plotting')
    #sim_time = time.time()

    ## Redshifts
    # Collapsar redshift pdf [number / yr] all-sky
    redshift_pdf_coll, N_coll = source_rate_density(redshifts,
        rho0=ps[par_dict["coll rho0"]],
        z_star=ps[par_dict["coll z*"]],
        n1=ps[par_dict["coll z1"]],
        n2=ps[par_dict["coll z2"]],
        vol_arr=vol_arr, plot=plot_func)
    # Merger redshift pdf [number / yr] all-sky
    redshift_pdf_merg, N_merg = source_rate_density(redshifts,
        rho0=ps[par_dict["merg rho0"]],
        z_star=ps[par_dict["merg z*"]],
        n1=ps[par_dict["merg z1"]],
        n2=ps[par_dict["merg z2"]],
        vol_arr =vol_arr, plot=plot_func)

    # Correct for 11.5 years of GBM data
    # Apparently, some distributions make N=inf,
    # so we have to ensure for that
    try:
        N_coll = np.int(N_coll*11.5)
        N_merg = np.int(N_merg*11.5)
    except:
        return -np.inf

    # Draw random redshifts from collapsar pdf
    redshift_sample_coll, _ = sample_distribution(redshifts,
        redshift_pdf_coll, num_draw=N_coll, xlabel='Redshift', ylog=False,
        plot=plot_func)
    # Draw random redshifts from merger pdf
    redshift_sample_merg, _ = sample_distribution(redshifts,
        redshift_pdf_merg, num_draw=N_merg, xlabel='Redshift', ylog=False,
        plot=plot_func)

    # Plot randomly drawn samples and pdf to ensure sim is done correctly
    if plot_func is not False:
        plot_samples(redshift_sample_coll, redshifts, redshift_pdf_coll)
        plot_samples(redshift_sample_merg, redshifts, redshift_pdf_merg)

    # Get kcorrection for each redshift
    coll_kc = kc[np.searchsorted(redshifts, redshift_sample_coll)]
    merg_kc = kc[np.searchsorted(redshifts, redshift_sample_merg)]
    # Get luminosity distance for each redshift
    coll_dl = dl[np.searchsorted(redshifts, redshift_sample_coll)]
    merg_dl = dl[np.searchsorted(redshifts, redshift_sample_merg)]


    ## Luminosity
    # Draw Collapsar GRB Rest-frame luminosities
    lum_pdf_coll, const_coll = get_luminosity(luminosities,
        lstar=1E52,
        alpha=ps[par_dict["coll alpha"]],
        beta=ps[par_dict["coll beta"]],
        lmin=luminosities[0], lmax=luminosities[-1], plot=plot_func)
    # Draw merger grb rest-frame luminosities
    lum_pdf_merg, const_merg = get_luminosity(luminosities,
        lstar=1E52,
        alpha=ps[par_dict["merg alpha"]],
        beta=ps[par_dict["merg beta"]],
        lmin=luminosities[0], lmax=luminosities[-1], plot=plot_func)

    # Draw random collapsar luminosities from pdf
    lum_sample_coll, lum_cdf_coll = sample_distribution(luminosities,
        lum_pdf_coll, num_draw=N_coll,
        xlabel='Peak Luminosity (1-10,000 keV)[ergs/s]',
        xlog=True, ylog=False, plot=plot_func)
    # Draw random merger luminosities from pdf
    lum_sample_merg, lum_cdf_merg = sample_distribution(luminosities,
        lum_pdf_merg, num_draw=N_merg,
        xlabel='Peak Luminosity (1-10,000 keV)[ergs/s]',
        xlog=True, ylog=False, plot=plot_func)
    '''
    I think this is ok, but I'm still super confused about the log(L)
    and normalization.
    # Plot randomly drawn samples and pdf to ensure sim is done correctly
    #if plot_func is not False:
    plot_samples(lum_sample_coll, luminosities, lum_pdf_coll, xlog=True,
        ylog=True)
    plot_samples(lum_sample_merg, luminosities, lum_pdf_merg, xlog=True,
        ylog=True)
    '''
    ## Peak flux
    # Get model peak flux for simulated collapsar GRBs
    coll_pf_all, coll_pf = Peak_flux(L=lum_sample_coll,
        z=redshift_sample_coll, kcorr=coll_kc, dl=coll_dl,
        emin=detector.getfloat('nai_emin'),
        emax=detector.getfloat('nai_emax'), plotting=plot_GRB)
    # Get model peak flux for simulated collapsar GRBs
    merg_pf_all, merg_pf = Peak_flux(L=lum_sample_merg,
        z=redshift_sample_merg, kcorr=merg_kc, dl=merg_dl,
        emin=detector.getfloat('nai_emin'),
        emax=detector.getfloat('nai_emax'), plotting=plot_GRB)
    
    # Combine collapsar and merger model counts
    model, data = combine_data(coll_model=coll_pf, merg_model=merg_pf,
        data=obs_pf, coll_all=coll_pf_all, merg_all=merg_pf_all,
        coll_model_label='Collapsar Model', merg_model_label='Merger Model',
        data_label='GBM Data', show_plot=plot_GRB)
    
    # Calculate the likelihood for this model (i.e., these parameters)
    try:
        llr = log_likelihood(model_counts=model, data_counts=data)
        
        # set uniform priors for right now
        ln_prior = 0
        for n in range(num_param):
            keyword = [x for x,y in par_dict.items() if y==n]
            ln_prior += np.log(prior_dist(ps[n], keyword[0]))
        return llr+ln_prior
        
    except:
        return -np.inf
