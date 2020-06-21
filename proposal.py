import numpy as np
from prior import *

"""
Pick new parameter
"""
def proposal_distribution(parameter, x, tuning=1):
    n = np.random.normal(x, tuning)
    if inbounds(n, parameter) is not True:
        # keep parameter at current value
        return -np.inf
    else:
        # return new value
        return n
'''
    if parameter == 'coll z1':
        # Low redshift index for Collapsars
        n = np.random.normal(x, tuning)
        if inbounds(n, parameter) is not True:
            # keep parameter at current value
            return -np.inf
        else:
            # return new value
            return n

    elif parameter == 'coll z2':
        # High redshift index for Collapsars
        n = np.random.normal(x, tuning)
        if inbounds(n, parameter) is not True:
            return -np.inf
        else:
            return n

    elif parameter == 'coll z*':
        # Collapsar Redshift peak
        n = np.random.normal(x, tuning*2)
        if inbounds(n, parameter) is not True:
            return -np.inf
        else:
            return  n

    elif parameter == 'coll rho0':
        # Local redshift rate for collapsars
        n = np.random.normal(x, tuning*2)
        if inbounds(n, parameter) is not True:
            return -np.inf
        else:
            return n

    elif parameter == 'merg z1':
        # Low redshift rate for mergers
        n = np.random.normal(x, tuning)
        if inbounds(n, parameter) is not True:
            return -np.inf
        else:
            return n

    elif parameter == 'merg z2':
        # Low redshift rate for mergers
        n = np.random.normal(x, tuning)
        if inbounds(n, parameter) is not True:
            return -np.inf
        else:
            return n

    elif parameter == 'merg z*':
        # Collapsar Redshift peak
        n = np.random.normal(x, tuning*2)
        if inbounds(n, parameter) is not True:
            return -np.inf
        else:
            return  n

    elif parameter == 'merg rho0':
        # Local redshift rate for collapsars
        n = np.random.normal(x, tuning*2)
        if inbounds(n, parameter) is not True:
            return -np.inf
        else:
            return n

    elif parameter == 'coll alpha':
        # Low luminosity index for collapsars
        n = np.random.normal(x, tuning)
        if inbounds(n, parameter) is not True:
            return -np.inf
        else:
            return n

    elif parameter == 'coll beta':
        # High luminosity index for mergers
        n = np.random.normal(x, tuning)
        if inbounds(n, parameter) is not True:
            return -np.inf
        else:
            return n

    elif parameter == 'merg alpha':
        # Low luminosity index for mergers
        n = np.random.normal(x, tuning)
        if inbounds(n, parameter) is not True:
            return -np.inf
        else:
            return n

    elif parameter == 'merg beta':
        # High luminosity index for mergers
        n = np.random.normal(x, tuning)
        if inbounds(n, parameter) is not True:
            return -np.inf
        else:
            return n

    elif parameter == 'coll mu':
        # Mean for collapsar duration distribution
        n = np.random.normal(x, tuning)
        if inbounds(n, parameter) is not True:
            return -np.inf
        else:
            return n
                
    elif parameter == 'merg beta':
        # Deviation from mean of collapsar duration distribution
        n = np.random.normal(x, tuning)
        if inbounds(n, parameter) is not True:
            return -np.inf
        else:
            return n
            
    elif parameter == 'merg beta':
        # High luminosity index for mergers
        n = np.random.normal(x, tuning)
        if inbounds(n, parameter) is not True:
            return -np.inf
        else:
            return n
            
    elif parameter == 'merg beta':
        # High luminosity index for mergers
        n = np.random.normal(x, tuning)
        if inbounds(n, parameter) is not True:
            return -np.inf
        else:
            return n
                    
    else:
        print ('parameter not found')
'''
