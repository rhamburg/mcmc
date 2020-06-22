import numpy as np

# Paramter bounds
#logarithmic?
collz1 = [0,5]
collz2 = [-5,0]
collzpeak = [0,10]
collrho0 = [0.01,10]
mergz1 = [0,5]
mergz2 = [-5,0]
mergzpeak = [0,10]
mergrho0 = [0.01, 10]
coll_alpha = [-5, 0]
coll_beta = [-5, 0]
merg_alpha = [-5, 0]
merg_beta = [-5, 0]
coll_mu = [1e-2, 100]
coll_sigma = [0, 10]
merg_mu = [1e-3, 10]
merg_sigma = [0, 10]
#coll_lpeak = [1e49, 1e54]
#merg_lpeak = [1e49, 1e54]


def check(x, a, b):
    """
    If parameter is within prior bounds, then return True

    Parameters
    ------------
    x: random number; can be int or float
    a: lower bound of prior distribution
    b: upper bound of prior distribution

    Returns
    ----------
    True/False - determined by whether or not x is within [a, b]
    """
    if x > a and x <= b:
        return True
    else:
        return False


def probability(x, a=0, b=1, type='uniform'):
    if type == 'uniform':
        return 1/(b-a)
    elif type == 'log uniform':
        return 1/(x*np.log(b/a))
    else:
        print ('you must specify another probability distribution')


def inbounds(n, parameter):
    """
    Check if parameter is within prior

    Parameters
    ____________
    n: random number drawn for proposed MCMC step; integer
    parameter: string describing which parameter is being randomly changed

    Returns
    ____________
    True/False - decided in the check function by whether n is within the bounds of its prior distribution
    """
    # Collapsar Redshift
    if parameter == 'coll z1':
        return check(n, collz1[0], collz1[1])
    elif parameter == 'coll z2':
        return check(n, collz2[0], collz2[1])
    elif parameter == 'coll z*':
        return check(n, collzpeak[0], collzpeak[1])
    elif parameter == 'coll rho0':
        return check(n, collrho0[0], collrho0[1])
    # Merger Redshift
    elif parameter == 'merg z1':
        return check(n, mergz1[0], mergz1[1])
    elif parameter == 'merg z2':
        return check(n, mergz2[0], mergz2[1])
    elif parameter == 'merg z*':
        return check(n, mergzpeak[0], mergzpeak[1])
    elif parameter == 'merg rho0':
        return check(n, mergrho0[0], mergrho0[1])
    # Collapsar Luminosity
    elif parameter == 'coll alpha':
        return check(n, coll_alpha[0], coll_alpha[1])
    elif parameter == 'coll beta':
        return check(n, coll_beta[0], coll_beta[1])
    # Merger Luminosity
    elif parameter == 'merg alpha':
        return check(n, merg_alpha[0], merg_alpha[1])
    elif parameter == 'merg beta':
        return check(n, merg_beta[0], merg_beta[1])
    # Collapsar duration
    elif parameter == 'coll mu':
        return check(n, coll_mu[0], coll_mu[1])
    elif parameter == 'coll sigma':
        return (check(n, coll_sigma[0], coll_sigma[1]))
    # Merger duration
    elif parameter == 'merg mu':
        return (check(n, merg_mu[0], merg_mu[1]))
    elif parameter == 'merg sigma':
        return (check(n, merg_sigma[0], merg_mu[1]))
    else:
        print ('parameter was not found; inbounds')


def prior_dist(n, parameter, type='uniform'):
    """
    Here we obtain the prior probability of each parameter n

    Parameters
    -----------
    n: the randomly drawn variable; int or float
    parameter: string that specified which parameter is being investigated
    type: type of known distribution; string

    Returns
    -----------
    the probability of the random variable n, given the prior; float
    """
    # Collapsar redshift
    if parameter == 'coll z1':
        return probability(n, a=collz1[0], b=collz1[1])
    elif parameter == 'coll z2':
        return probability(n, a=collz2[0], b=collz2[1])
    elif parameter == 'coll z*':
        return probability(n, a=collzpeak[0], b=collzpeak[1])
    elif parameter == 'coll rho0':
        return probability(n, a=collrho0[0], b=collrho0[1], type='log uniform')
    # Merger redshift
    elif parameter == 'merg z1':
        return probability(n, a=mergz1[0], b=mergz1[1])
    elif parameter == 'merg z2':
        return probability(n, a=mergz2[0], b=mergz2[1])
    elif parameter == 'merg z*':
        return probability(n, a=mergzpeak[0], b=mergzpeak[1])
    elif parameter == 'merg rho0':
        return probability(n, a=mergrho0[0], b=mergrho0[1], type='log uniform')
    # Collapsar luminosity
    elif parameter == 'coll alpha':
        return probability(n, a=coll_alpha[0], b=coll_alpha[1])
    elif parameter == 'coll beta':
        return probability(n, a=coll_beta[0], b=coll_beta[1])
    # Merger luminosity
    elif parameter == 'merg alpha':
        return probability(n, a=merg_alpha[0], b=merg_alpha[1])
    elif parameter == 'merg beta':
        return probability(n, a=merg_beta[0], b=merg_beta[1])
    # Collapsar duration
    elif parameter == 'coll mu':
        return probability(n, a=coll_mu[0], b=coll_sigma[1], type='log uniform')
    elif parameter == 'coll sigma':
        return probability(n, a=coll_sigma[0], b=coll_sigma[1])
    # Merger duration
    elif parameter == 'merg mu':
        return probability(n, a=merg_mu[0], b=merg_mu[1], type='log uniform')
    elif parameter == 'merg sigma':
        return probability(n, a=merg_sigma[0], b=merg_sigma[1])
    else:
        print ('parameter was not found; prior_dist')
