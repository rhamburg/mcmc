import numpy as np

# Paramter bounds
collz1 = [0,5]
collz2 = [-5,0]
collzpeak = [0,10]
collrho0 = [0,10] #logarithmic?
mergz1 = [0,5]
mergz2 = [-5,0]
mergzpeak = [0,10]
mergrho0 = [0,10] #logarithmic?
coll_alpha = []
coll_beta = []
coll_lpeak = []
merg_alpha = []
merg_beta = []
merg_lpeak = []


def check(x, a, b):
    """ If parameter is within prior bounds, then return True"""
    if x > a and x < b:
        return True
    else:
        return False

def probability(x, a=0, b=1, type='uniform'):
    if type == 'uniform':
        return 1/(b-a)
    else:
        print ('you must specify another probability distribution')

def inbounds(n, parameter):
    """Check if parameter is within prior"""
    if parameter == 'coll z1':
        return check(n, collz1[0], collz1[1])
    elif parameter == 'coll z2':
        return check(n, collz2[0], collz2[1])
    elif parameter == 'coll z*':
        return check(n, collzpeak[0], collzpeak[1])
    elif parameter == 'coll rho0':
        return check(n, collrho0[0], collrho0[1])
    else:
        print ('parameter was not found')

def prior_dist(n, parameter, type='uniform'):
    if parameter == 'coll z1':
        return probability(n, a=collz1[0], b=collz1[1])
    elif parameter == 'coll z2':
        return probability(n, a=collz2[0], b=collz2[1])
    elif parameter == 'coll z*':
        return probability(n, a=collzpeak[0], b=collzpeak[1])
    elif parameter == 'coll rho0':
        return probability(n, a=collrho0[0], b=collrho0[1])
    else:
        print ('parameter was not found')
