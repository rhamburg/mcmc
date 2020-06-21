from astropy.io import ascii
from astropy.cosmology import FlatLambdaCDM
import numpy as np
from numpy import random
import matplotlib.pyplot as plt
import scipy
from scipy.constants import c
from scipy.integrate import quad, quad_explain
from math_functions import *
import matplotlib
from plot import *
import time


def read_data(type='all'):
    """
    # GBM's peak flux and duration (t90) distributions
    
    Parameters:
    ------------
    type: string of 'long', 'short', or 'all'
    
    Returns:
    ------------
    peak fluxes, durations lists
    """
    file = ascii.read("gbm_bursts_t90_flux.txt", data_start=1)
    if type == 'long':
        long_GRBs = file[np.where(file["col3"]>2.)]
        return long_GRBs["col5"], long_GRBs["col3"]
    elif type == 'short':
        short_GRBs = file[np.where(file["col3"]<2.)]
        return peak_flux_long["col5"], short_GRBs["col3"]
    else:
        return file["col5"], file["col3"]

"""
# Find cdf interval j such that the random number is in that interval
# Draw randomly from a distribution
Parameters:
-----------
xpdf: pdf of distribution you want to draw from
cdf: cdf of distribution you want to draw from
num_draw: num of random draws from the distribution
"""
def draw_samples(pdf, cdf, num_draw=1000):
    #random.seed(3)
    draws = random.random(num_draw)
    sample = [pdf[np.searchsorted(cdf, draws)]][0]
    return sample


def get_luminosity(lum, alpha=0.2, beta=1.4, lstar=1E52, lmin=1E49, lmax=1E52,
plot=False):
    """
    # Intrinsic luminosity distribution (number of GRBs/luminosity interval)

    Parameters:
    -----------
    lum: input luminosities
    alpha: low luminosity index
    beta: high luminosity index
    lstar: break of luminosity function
    lmin: lower luminosity limit from which we integrate luminosity function
    lmax: upper luminosity limit to which we integrate luminosity function

    Returns:
    pdf_normed: the normalized pdf of the luminosity function
    A: the integration constant
    """
    # Integrate the luminosity function to find the normalization constant
    sum = quad(luminosity_function, lmin, lmax, args=(lstar, alpha, beta, 1.))
    # integration constant
    A = 1 / sum[0]
    # Obtain luminosity pdf
    pdf_normed = [luminosity_function(l, phi_0=A, lum_star=lstar, \
        alpha=alpha, beta=beta) for l in lum]
    if plot is not False:
        plt.xscale('log')
        plt.yscale('log')
        plt.plot(lum, pdf_normed)
        plt.xlabel('Peak luminosity (1-10,000 keV) [ergs/s]')
        plt.ylabel('dN/dlogL')
        plt.grid()
        plt.show()
        plt.close()
    return pdf_normed, A


def sample_distribution(xs, pdf, num_draw=1000, plot=False, xlabel=None,
xlog=False, ylog=False):
    """
    # Sample luminosity and redshift distributions

    Parameters:
    -----------
    xs: input depedent values
    pdf: distribution pdf
    num_draw: number of GRBs to randomly draw from cdf
    plot: option to plot cdf
    xlabel: xlabel for plot
    xlog: option for plotting x-axis in logscale
    ylog: option for plotting y-axis in logscale

    Returns:
    -----------
    random_xs: the randomly drawn xs's from the pdf
    cdf_normed: the normalized cdf of the pdf
    """
    cdf = np.cumsum(pdf)
    cdf_normed = [c/cdf[-1] for c in cdf]
    random_xs = draw_samples(xs, cdf_normed, num_draw=num_draw)
    if plot is not False:
        plot_cdf(xs, cdf_normed, xlabel, xlog, ylog)
    return random_xs, cdf_normed


def source_rate_density(redshifts, rho0=1., z_star=1., n1=1., n2=1.,
vol_arr=None, plot=False):
    """
    # Calculates 'observed' GRB redshift distribution

    Parameters:
    -----------
    redshifts: array of redshifts to calculate rate at
    rho0: rate at z=0; local rate density ρ0 (in Gpc− 3 yr− 1)
    z_star: redshift at which the rate peaks and turns over
    n1: the low redshift index (positive)
    n2: the high redshift index (negative)
    vol_arr: array of precalculated comoving volume elements for each redshift
        in the redshifts array
    plot: Decision to plot or not

    Returns:
    -----------
    rate_pdf: the normalized rate of GRBs / yr / dz
    N_tot: the total number of GRBs / yr in the entire universe
    """
    # Number of GRBs in the universe per year (no jet angle included)
    N_tot = quad(obs_grb_rate, redshifts[0], redshifts[-1],
        args=(z_star, n1, n2, rho0, None, False))
    # Observed rate of GRBs / year / dz
    rate = np.array([obs_grb_rate(redshifts[i], z_star=z_star, n1=n1, n2=n2, \
        r_0=rho0, vol_arr=vol_arr[i], precalc=True) \
        for i in range(len(redshifts))])
    # Normalize rate by total number of GRBs / yr
    rate_pdf = rate / N_tot[0]
    # Plot 'observed' redshift distribution
    if plot is not False:
        int_rate = [int_grb_rate(z, z_star=z_star, n1=n1, n2=n2, r_0=rho0) \
            for z in redshifts]
        plot_rate(redshifts, int_rate, title='Intrinsic GRB Rate',
            ylabel=r'$R_{GRB;int}$ [$Gpc^{-1} yr^{-1}$]')
        plot_rate(redshifts, rate, ylabel=r'$R_{GRB;obs}$ [$dz^{-1} yr^{-1}$]',
            title='Obs GRB Rate')
    return rate_pdf, N_tot[0]


"""
Calculate peak flux
"""
def Peak_flux(L=None,z=None, threshold=1., kcorr=None, emin=10, emax=1000,
dl=None, plotting=False):
    # Peak photon flux
    pf = (1+z) * (L / (4 * np.pi * dl**2)) * kcorr
    corr_pf = np.array(pf)
    # Correct for GBM FOV
    if len(corr_pf) > 0:
        num_fov = int(len(corr_pf) * 0.67 * 0.85)
        idx_fov = random.random_integers(0, len(corr_pf)-1, num_fov)
        corr_pf = corr_pf[idx_fov]
        L = L[idx_fov]
        z = z[idx_fov]
        # Take peak fluxes above detection limit
        detected = det_thresh(corr_pf, kind='prob')
        final_pf = corr_pf[detected]
        final_l = L[detected]
        final_z = z[detected]
        # Plot stuff
        if plotting is not False:
            plot_peak_flux_color(final_l, final_z, final_pf)
            plot_obs_lum(final_l)
            plot_obs_rate(final_z)
        return corr_pf, final_pf
    else:
        return [], []


"""
Detector threshold
"""
def det_thresh(pf, threshold=0.5, kind='cutoff', plot=True):
    if kind == 'prob':
        # Obtain detection probability for give peak flux
        det_prob = logistic(pf, s=0.1, w=5)
        # If random sample from binomial probability
        # yields value greater than detection probability
        # then accept, else reject
        binary = np.random.binomial(1, det_prob)
        return (binary == 1)
    elif kind == 'cutoff':
        return np.where(pf[:]>threshold)


def log_likelihood(model_counts, data_counts):
    """
    # Calculate the log likelihood
    # When calculating Poisson probability, we ignore the factorial term because
    # it's model independent

    Parameters:
    -----------
    model_counts: normalized or unnormalized counts for model
    data_counts: normalized or unnormalized counts for data

    Returns:
    -----------
    LLR - float; log likelihood calculated from model counts and data counts
    """
    # Replace 0's in model rate with something very small
    model_counts = model_counts.astype(float)
    model_counts[model_counts==0] = 1e-20
    LLR = np.sum(-model_counts + data_counts*np.log(model_counts))
    return LLR


def combine_data(coll_model=None, merg_model=None, data=None, coll_all=None,
merg_all=None, show_plot=False, coll_model_label=None, merg_model_label=None,
data_label=None):
    """
    Combine data to prepare for LLR calculation

    Parameters:
    ------------
    coll_model - array or list of "detected" peak fluxes from model GRBs from collapsars
    merg_model - array or list of "detected" peak fluxes from model GRBs from mergers
    data - array of peak fluxes from real data
    coll_all - array or list of all (detected and undetected) model peak fluxes from collapsars
    merg_all - array of list of all (detected and undetected) model peak fluxes from mergers
    coll_model_label -
    merg_model_label -
    data_label -

    Returns:
    ------------
    model_counts - list of GRBs in each collapsar+merger peak flux histogram bin
    data_counts - list of GRBs in each data peak flux histogram bin
    """
    tot_model = np.concatenate([coll_model, merg_model])
    bins=np.logspace(-1, 3, 70)
    model_counts, model_bin_edges = np.histogram(tot_model, bins=bins)
    data_counts, data_bin_edges = np.histogram(data, bins=bins)
    if show_plot is not False:
        plot_peak_flux(coll_all, merg_all, coll_model, merg_model, data,
            coll_model_label, merg_model_label, tot_model)
    return model_counts, data_counts


"""
Metro-hastings decision algorithm
"""
def metro_hastings(current_post=None, proposed_post=None, current_param=None,
p_accept=None, total_accept=None, proposed_param=None):
    rand_num = np.random.uniform(0,1)
    H = np.minimum(1, proposed_post / current_post)
    if H >= rand_num and np.isinf(proposed_post)==False:
        # accept candidate
        # update acceptance ratios
        #print('accept')
        total_accept += 1
        p_accept += 1
        return proposed_param, proposed_post, p_accept, total_accept
    else:
        # reject candidate
        # return to previous parameter value and posterior
        return current_param, current_post, p_accept, total_accept

"""
Metro-hastings decision algorithm
"""
def metro_hastings_test(current_val, proposed_val, current_post, proposed_post):
    rand_num = np.random.uniform(0,1)
    H = np.minimum(1., np.exp(proposed_post - current_post))
    #print(H, rand_num)
    if H >= rand_num:
        # accept candidate
        print ('accept')
        #total_accept += 1
        # update acceptance ratios
        #p_accept += 1
        #print (proposed_val, proposed_post)
        return proposed_val, proposed_post
    else:
        # reject candidate
        print ('reject')
        # return to previous parameter value and posterior
        #print (current_val, current_post)
        return current_val, current_post

'''
# Intrinsic durations distributions
def intrinsic_durations(N, mint=-1, maxt=3, mmint=-3, mmaxt=2, coll_mid=2, merg_mid=0.1, coll_sig=1, merg_sig=1, plot=False):
    # Collapsars
    collapsar_t90s = np.logspace(mint,maxt,N)
    coll_mu = np.log(coll_mid) + coll_sig**2
    csum = quad(log_normal, 10**mint, 10**3, args=(coll_mu, coll_sig))
    coll_t90_freq = log_normal(collapsar_t90s, mu=coll_mu, sigma=coll_sig, A=1./csum[0])
    # Mergers
    merger_t90s = np.logspace(mmint,mmaxt, N)
    merg_mu = np.log(merg_mid) + merg_sig**2
    msum = quad(log_normal, 10**mmint, 10**mmaxt, args=(merg_mu, merg_sig))
    # Add scaling factor of 1/20 just to normalize distributions to themselves
    merger_t90_freq = log_normal(merger_t90s, mu=merg_mu, sigma=merg_sig, A=1./(20*msum[0]))
    # Sample distributions
    ct90 = random.lognormal(coll_mu, coll_sig, N*10)
    mt90 = random.lognormal(merg_mu, merg_sig, N*10)
    ctbins = np.logspace(mint,maxt,50)
    mtbins = np.logspace(mmint,mmaxt,50)
    if plot is not False:
        matplotlib.rc('xtick', labelsize=14)
        matplotlib.rc('ytick', labelsize=14)
        plt.plot(collapsar_t90s, coll_t90_freq, linewidth=2, label='Collapsars')
        #plt.hist(ct90, bins=ctbins, density=True, color='#2ca02c', label='Collapsar random')
        plt.plot(merger_t90s, merger_t90_freq, linewidth=2, label='Mergers')
        #plt.hist(mt90, bins=mtbins, density=True, color='#d62728', label='Mergers random')
        plt.xscale('log')
        plt.xlabel(r'Duration [s]',fontsize=16)
        plt.ylabel('PDF',fontsize=16)
        plt.legend(fontsize=16)
        #plt.show()
        plt.tight_layout()
        plt.savefig('new_durations',dpi=300)
        plt.close()
    return collapsar_t90s, coll_t90_freq, ct90, mt90

# Simulate a light curve
def lightcurve(times, duration, rand_lum, limit=0.01, plot=False):
    # Norris lightcurve
    pulse = [norris(t, tau1=39, tau2=4) for t in times]
    # Adjust the width to fit the random duration of the burst
    a = duration / 50.
    rest_times = [a*t for t in times]
    # Find the start (rise to 1% above background)
    for i in range(len(pulse)):
        if pulse[i] > limit:
            start_time = rest_times[i]
            index = i
            break
    # Find the end (return to 1% above background)
    for i in range(index+1, len(pulse)):
        if pulse[i] < limit:
            end_time = rest_times[i]
            break
    #print ('start', start_time)
    #print ('end', end_time)
    # Adjust amplitude to duration
    pulse = [p*rand_lum for p in pulse]
    if plot is not False:
        # Plot GRB lightcurve
        plt.axhline(y=0.0, color='black', linestyle='-')
        plt.axvline(x=start_time, color='r', linestyle='--')
        plt.axvline(x=end_time, color='r', linestyle='--')
        plt.plot(rest_times, pulse)
        plt.ylabel('Luminosity [ergs/s]')
        plt.xlabel('Time [s]')
        plt.title('Simulated GRB')
        plt.show()
        plt.close()
    return (pulse, rest_times, start_time, end_time)

# Find the observed duration
def find_obs_dur(times, pf, limit=0.8, plot=False):
    end = None
    for f in range(len(pf)):
        if pf[f] > limit:
            start = times[f]
            index = f
            break
    for f in range(index+1, len(pf)):
        if pf[f] < limit:
            end = times[f]
            new_inx = f
            break
    if end is None:
        end = times[-1]
        new_inx = len(times)-1

    # Fluence
    new_pf = pf[index:new_inx]
    new_t = times[index:new_inx]
    fluence = 0
    for l in range(len(new_pf)-1):
        fluence += new_pf[l] * (times[l+1] - times[l])
    t90 = end - start
    return t90, fluence
'''
