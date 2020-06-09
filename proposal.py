import numpy as np
from prior import *

"""
Pick new parameter
"""
def proposal_distribution(parameter, x, tuning=1):
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
