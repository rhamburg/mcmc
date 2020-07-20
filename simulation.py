import numpy as np
from numpy import random
from prior import prior_dist
from source_functions import *
from plot import plot_samples
import tempfile
import time

def Simulation(ps, par_dict, num_param, redshifts=None, luminosities=None,
durations=None, obs_pf=None, obs_t90=None, detector=None, options=None,
vol_arr=None, kc=None, dl=None, plot_GRB=None, prior=False, sim=False, dsim=False,
file=None):

    # Option for prior testing
    if prior is not False:
        ln_prior = 0
        for n in range(num_param):
            keyword = [x for x,y in par_dict.items() if y==n]
            ln_prior += np.log(prior_dist(ps[n], keyword[0]))
        return ln_prior

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
    print (N_coll, N_merg)
    
    # Draw random redshifts from collapsar pdf
    redshift_sample_coll = sample_distribution(redshifts,
        redshift_pdf_coll, xlabel='Redshift', ylog=False,
        num_draw=N_coll, plot=plot_func)
    # Draw random redshifts from merger pdf
    redshift_sample_merg = sample_distribution(redshifts,
        redshift_pdf_merg, xlabel='Redshift', ylog=False,
        num_draw=N_merg, plot=plot_func)

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
    lum_sample_coll = sample_distribution(luminosities,
        lum_pdf_coll, xlabel='Peak Luminosity (1-10,000 keV)[ergs/s]',
        num_draw=N_coll, xlog=True, ylog=False, plot=plot_func)
    # Draw random merger luminosities from pdf
    lum_sample_merg = sample_distribution(luminosities,
        lum_pdf_merg, xlabel='Peak Luminosity (1-10,000 keV)[ergs/s]',
        num_draw=N_merg, xlog=True, ylog=False, plot=plot_func)
    
    # Plot randomly drawn samples and pdf to ensure sim is done correctly
    if plot_func is not False:
        bins = np.logspace(50, 54, 60)
        plot_samples(lum_sample_coll, luminosities, lum_pdf_coll, bins=bins, xlog=True, ylog=True)
        plot_samples(lum_sample_merg, luminosities, lum_pdf_merg, bins=bins, xlog=True, ylog=True)
        
    
    ## Peak flux
    # Get model peak flux for simulated collapsar GRBs
    coll_pf_all, coll_pf, coll_z = Peak_flux(L=lum_sample_coll,
        z=redshift_sample_coll, kcorr=coll_kc, dl=coll_dl, plotting=plot_GRB,
        sim=sim, dsim=dsim, title='Collapsars')
    # Get model peak flux for simulated collapsar GRBs
    merg_pf_all, merg_pf, merg_z = Peak_flux(L=lum_sample_merg,
        z=redshift_sample_merg, kcorr=merg_kc, dl=merg_dl, plotting=plot_GRB,
        sim=sim, dsim=dsim, title='Mergers')

    '''
    ## Duration
    # Collapsar duration pdf
    coll_dur_pdf = intrinsic_duration(durations, mu=ps[par_dict["coll mu"]], sigma=ps[par_dict["coll sigma"]], plot=plot_func)
    dur_sample_coll = sample_distribution(durations, coll_dur_pdf, num_draw=len(coll_pf), plot=plot_func)

    # Merger duration pdf
    merg_dur_pdf = intrinsic_duration(durations, mu=ps[par_dict["merg mu"]], sigma=ps[par_dict["merg sigma"]], plot=plot_func)
    dur_sample_merg = sample_distribution(durations, merg_dur_pdf, num_draw=len(merg_pf), plot=plot_func)
    
    # Plot randomly drawn samples and pdf to ensure sim is done correctly
    #if plot_func is not False:
    bins = np.logspace(-2, 3, 60)
    plot_samples(dur_sample_coll, durations, coll_dur_pdf, bins=bins, xlog=True)
    plot_samples(dur_sample_merg, durations, merg_dur_pdf, bins=bins, xlog=True)
    
    
    #Also need to adjust duration energy range to match t90 somehow...
    
    # Observed duration is longer than source duration
    dur_sample_coll *= (coll_z + 1)
    dur_sample_merg *= (merg_z + 1)
    '''
    # Combine collapsar and merger model counts
    pf_model, pf_data = combine_data(coll_model=coll_pf,merg_model=merg_pf,
        data=obs_pf, coll_all=coll_pf_all, merg_all=merg_pf_all,
        coll_model_label='Collapsar Model', merg_model_label='Merger Model',
        data_label='GBM Data', show_plot=plot_GRB, sim=sim)
    
    # Combine collapsar and merger model durations
    #dur_model, dur_data = combine_data(coll_model=dur_sample_coll, merg_model=dur_sample_merg, data=obs_t90, coll_all=dur_sample_coll, merg_all=dur_sample_merg, show_dur_plot=True)#plot_GRB)

    # Save peak flux data
    if sim is not False and file is not None:
        print ('Saving simulation data: ../'+file+'.npy')
        tot_model = np.concatenate([coll_pf, merg_pf])
        np.save('../'+file+'.npy', tot_model)
        

    # Calculate the likelihood for this model (i.e., these parameters)
    try:
        pf_llr = log_likelihood(model_counts=pf_model, data_counts=pf_data)
        dur_llr = 0#log_likelihood(model_counts=dur_model, data_counts=dur_data)
      
        # set uniform priors for right now
        ln_prior = 0
        for n in range(num_param):
            keyword = [x for x,y in par_dict.items() if y==n]
            ln_prior += np.log(prior_dist(ps[n], keyword[0]))
        return pf_llr+dur_llr+ln_prior

    except:
        return -np.inf
