import math
import numpy as np
from numpy import random
from scipy.constants import c
from astropy.cosmology import FlatLambdaCDM
from scipy.integrate import quad
from scipy.special import erf

# Estimate normalized autocorrelation function of the stochastic process
# that generated the given dataset
# Parameters:
#   x - dataset array
#   lags - lag array
''' Need to write as FFT'''
def autocorrelation(x, lags):
    xmean = np.mean(x)
    N = len(x)
    rho = []
    for k in lags:
        num = (1./(N-k)) * np.sum((x[:N-k] - xmean) * (x[k:] - xmean))
        denom = (1./N) * np.sum((x - xmean)**2)
        rho.append(num/denom)
    return np.array(rho)

# Estimate the autocorrelation time,
# given an estimate of the autocorrelation function
# rho - array of autocorrelation values for lags
# M - ; should be << N
def autocorrelation_time(rho, M):
    return (1 + 2 * np.sum(rho[:M]))

# Symmetric logistic curve for a log-x domain
# x - data
# w - stretch/compression factor; larger w means higher compression
# s - shift factor; larger s means shifting to the left
def logistic(x, s=1, w=1):
    return (s*(x**w)) / (1 + s*(x**w))

# Binomial distribution
# p - Probability
# n - total trials
# k - number of successes
def binomial(p, n, k):
    nchoosek = math.factorial(n) / (math.factorial(k) * math.factorial(n-k))
    return nchoosek * (p**k) * ((1-p)**(n-k))

# Gaussian distribution
def gaussian(x, mu=0, sigma=1):
    return 1./ (sigma * np.sqrt(2*np.pi)) * np.exp( - (x-mu)**2 / (2 * sigma**2))

# Log normal function
def log_normal(x, mu=0, sigma=1, A=1):
    return A * np.exp(-(np.log(x) - mu)**2 / (2 * sigma**2))/ (x * sigma * np.sqrt(2 * np.pi))

# Log normal cdf
def log_normal_cdf(x, mu=0, sigma=1.):
    z = (np.log(x) - mu) / (np.sqrt(2) * sigma)
    return 0.5 + (0.5 * erf(z))

# Norris pulse shape
def norris(t, tau1=1, tau2=1, A=1):
    mu = (tau1/tau2)**(1/2)
    lamb = np.exp(2*mu)
    intensity = A * lamb * np.exp((-tau1/t) - (t/tau2))
    return intensity

def band(E, alpha=-1, beta=-2.3, E_0=511, N_0=1., energy_units=False):
    """
    INFO
    The Band function [photons cm^-2 keV^-1]
    Input energy should be in keV
    alpha - low energy index; beta - high energy index
    E_p = (2+alpha)E_0 - it's the peak energy in units of keV in the observed
        nuFnu spectrum
     N_0 is normalization factor at 100 keV [photons cm^-2 keV^-1]

    FOR TEST CASES,
    First integral of band, type in WolframAlpha:
       integral from 10 to 511  (x/100)^(-1) * e^(-x/511) dx
    For second integral, type
        integral from 511 to 1000  (511/100) * e^(-1) * (x/100)^(-2) dx
    and compare results with what numpy quad gives
    """
    E_b = (alpha - beta)*E_0 # 511 with initial parameters
    if E <= E_b:
        N = N_0 * ((E/100.)**alpha) * np.exp(-E/E_0)
    else:
        N = N_0 * ((E_b/100.)**(alpha-beta)) * np.exp(beta-alpha) * ((E/100)**beta)
    if energy_units is not False:
        return E*N
    else:
        return N

# Error function
def integral(t):
    return np.exp(-t**2)
def error_func(z):
    return 2./np.sqrt(np.pi)*(quad(integral, 0, z))[0]

# Calculate the k-correction
def k_correction(z, emin=10, emax=1000):
    num = quad(band, emin, emax)
    denom = quad(band, emin*(1+z), emax*(1+z))
    c_k = num[0] / denom[0]
    return c_k

# Luminosity function
def luminosity_function(l, lum_star=1E52, alpha=0.2, beta=1.4, phi_0=1.):
    if l < lum_star:
        index = alpha
    else:
        index = beta
    phi = phi_0 * (l/lum_star)**index
    return phi

# Intrinsic GRB rate function
# Should follow SFR?
# units of [M (Mpc)^-3 yr^-1]
def int_grb_rate(z, z_star=3., n1=2., n2=-2., r_0=1.):
    if z <= z_star:
        r = r_0 * (1+z)**n1
    else:
        r = r_0 * (1+z_star)**(n1-n2) * (1+z)**n2
    return r

# Observed GRB rate function [ M Gpc^-3 yr^-1]
# local rate at z=0: r_0 [Gpc^−3 yr^−1]
#   integral of rgrb from 0-3 = 21
#   integral of rgrb from 3-10 = 64.742...
#   integral of 1/(1+z) = 2.3978...
#   integral of rgrb / (1+z) = 17.68....
def obs_grb_rate(z, z_star=3., n1=2., n2=-2., r_0=1., vol_arr=None, precalc=True):
    rgrb = int_grb_rate(z, z_star=z_star, n1=n1, n2=n2, r_0=r_0)
    if precalc is not True:
        return diff_comoving_volume(z) * rgrb / (1+z)
    else:
        return vol_arr * rgrb / (1+z)

# Comoving volume element
# units of number per solid angle (dΩ), per redshift interval (dz), and per time
# interval in the observed frame (dtobs)
# dV / dz
def diff_comoving_volume(z):
    cosmo = FlatLambdaCDM(H0=70, Om0=0.3, Tcmb0=2.725)
    H0 = cosmo.H(0).value  # km / s / Mpc
    dl = cosmo.luminosity_distance(z).value
    lightspeed = c * 0.001 # convert from m/s to km/s
    ez = np.sqrt( cosmo.Om(0)*((1+z)**3) + (1-cosmo.Om(0)) )
    dv = 4*np.pi * (lightspeed/H0) * (dl**2) / ((1+z)**2) / ez
    dv = dv * ((0.001)**3) # Mpc^3 to Gpc^3
    # Test case
    #dv_new = cosmo.differential_comoving_volume(3)
    #print (4 * np.pi * dv_new.value)
    return dv
