import numpy as np

# Parameter dictionary
parameter_bounds = {
"coll z1":    [0, 5, 'uniform'],
"coll z2":    [-5, 0, 'uniform'],
"coll z*":    [1e-5, 6, 'log-uniform'],
"coll rho0":  [0.01,10, 'log-uniform'],
"merg z1":    [0, 5, 'uniform'],
"merg z2":    [-5, 0, 'uniform'],
"merg z*":    [1e-5, 6, 'log-uniform'],
"merg rho0":  [0.01,10, 'log-uniform'],
"coll alpha": [-5, 0, 'uniform'],
"coll beta":  [-5, 0, 'uniform'],
"merg alpha": [-5, 0, 'uniform'],
"merg beta":  [-5, 0, 'uniform'],
"coll mu":    [1e-3,100, 'log-uniform'],
"coll sigma": [1e-3, 10, 'log-uniform'],
"merg mu":    [1e-3, 10, 'log-uniform'],
"merg sigma": [1e-3, 10, 'log-uniform']
}
#coll_lpeak = [1e49, 1e54]
#merg_lpeak = [1e49, 1e54]


def proposal_distribution(parameter, x, tuning=1):
    """
    Pick new parameter
    
    Parameters:
    -------------
    parameter: string for the parameter to pick
    x: float, current value for the parameter
    tuning: float, width/variance of the proposal distribution
    
    Returns:
    -------------
    n: float, either the new value of the parameter of -inf if the selection was
        outside parameter's prior bounds
    """
    n = np.random.normal(x, tuning)
    if inbounds(n, parameter) is not True:
        # keep parameter at current value
        return -np.inf
    else:
        # return new value
        return n
        
        
def check(x, a=None, b=None):
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
    elif type == 'log uniform' or type=='log-uniform':
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
    try:
        return check(n, a=parameter_bounds.get(parameter)[0], b=parameter_bounds.get(parameter)[1])
    except:
        print ('parameter was not found in inbounds()')
        

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
    try:
        return probability(n, a=parameter_bounds.get(parameter)[0], b=parameter_bounds.get(parameter)[1], type=parameter_bounds.get(parameter)[2])
    except:
        print ('parameter was not found in prior_dist()')
        
