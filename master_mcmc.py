from argparse import ArgumentParser
from configparser import ConfigParser
from astropy.cosmology import FlatLambdaCDM
from math_functions import band
import numpy as np
from prior import prior_dist
from proposal import proposal_distribution
from scipy.integrate import quad
from simulation import Simulation
from source_functions import diff_comoving_volume, metro_hastings, metro_hastings_test, read_data
import time

# Read arguments
parser = ArgumentParser(prog='Thesis Program', description='Fitting GBM Peak Flux Distributions')
parser.add_argument('-f', '--filename', help='Filename to save results under')
parser.add_argument('-i', '--iter', help='Number of MCMC iterations to make', type=int)
parser.add_argument('-p', '--plotGRB', help='Boolean for plotting peak flux info')
parser.add_argument('-n', '--num_param', help='Number of parmeters to search')
parser.add_argument('-c', '--config', help='Path to config file to initialize parameters')
parser.add_argument('-P', '--prior', help='Boolean to determine if just sampling priors and not posterior')
args = parser.parse_args()


if args.iter is None:
    args.iter = 0
if args.plotGRB is not None:
    plot = args.plotGRB
else:
    plot = None
if args.num_param is None:
    raise Exception('Need to specifiy number of parameters to search')
else:
    num_param = int(args.num_param)
if args.filename is not None:
    file = 'results/'+args.filename+'.npy'
if args.config is not None:
    config_file = args.config
else:
    config_file = 'parameters.ini'
if args.prior is None:
    args.prior = False
else:
    print ('TESTING PRIOR DISTRIBUTIONS...')

config = ConfigParser()
config.read(config_file)
iterations = config['iterations']
detector = config['detector']
luminosity = config['luminosity']
redshift = config['redshift']
options = config['options']


print ('\nGRB Peak Flux Simulator\n')


# Read GBM data [peak flux 10-1000 keV]
print ('Reading data file...')
obs_pf = read_data()

# Luminosity info
NSIM = iterations.getint('NSIM')
min_lum = luminosity.getfloat('coll_min_lum')
max_lum = luminosity.getfloat('coll_max_lum')
luminosities = np.logspace(min_lum, max_lum, NSIM)

# Redshift info
zmin = redshift.getfloat('min_z')
zmax = redshift.getfloat('max_z')
redshifts = np.linspace(zmin, zmax, NSIM)

# Calculate differential comoving volume
print ('Calculating differential comoving volume...')
v_comov = diff_comoving_volume(redshifts)

# Calculate luminosity distance
print ('Calculating luminosity distance...')
cosmo = FlatLambdaCDM(H0=70, Om0=0.3, Tcmb0=2.725)
dl = cosmo.luminosity_distance(redshifts)
dl = dl.value*3.086e+24 # convert Mpc to cm

print ('Calculating detector correction and k-correction...')
# Spectral info
emin=detector.getint('nai_emin')
emax=detector.getint('nai_emax')
E_0 = detector.getint('E_0')
alpha = detector.getfloat('alpha')
beta = detector.getfloat('beta')
# Detector correction
k_num = quad(band, emin, emax, args=(alpha, beta, E_0, 1., False))[0]
# K-correction
k_denom = np.array([quad(band, 1/(1+z), 10000/(1+z),
    args=(alpha, beta, E_0, 1., True))[0] for z in redshifts])
# (converted from spectral photon flux to spectral energy flux)
k_denom *= 1.6e-9
kcorr = k_num/k_denom



print ('\nStarting mcmc...')
print ('FIX LUMINOSITY??')
sim_time = time.time()

# Initialize parameter space
parameter_dict = {
    "coll z1": 0,
    "coll z2": 1,
    "coll z*": 2,
    "coll rho0": 3,
    "merg z1": 4,
    "merg z2": 5,
    "merg z*": 6,
    "merg rho0": 7,
    "coll alpha": 8,
    "coll beta":  9,
    "merg alpha": 10,
    "merg beta": 11,
    "posterior": 12,
    "acceptance num": 13,
    "effective step": 14,
    }

# some easy index definitions used in MCMC loop
post_idx = parameter_dict["posterior"]
p_accept_idx = parameter_dict["acceptance num"]
eff_step_idx = parameter_dict["effective step"]

# Set up parameter space to record MCMC inputs
rows, cols = (args.iter+1, len(parameter_dict))
parameter_space = np.array([[0.0]*cols]*rows)
parameter_space[0] = [redshift.getfloat('coll_n1'),
                        redshift.getfloat('coll_n2'),
                        redshift.getfloat('coll_zstar'),
                        redshift.getfloat('coll_rho0'),
                        redshift.getfloat('merg_n1'),
                        redshift.getfloat('merg_n2'),
                        redshift.getfloat('merg_zstar'),
                        redshift.getfloat('merg_rho0'),
                        luminosity.getfloat('coll_alpha'),
                        luminosity.getfloat('coll_beta'),
                        luminosity.getfloat('merg_alpha'),
                        luminosity.getfloat('merg_beta'), 1., 0., 0]

'''
# Loop through each tuning parameter
burn_out = False
if burn_out is not False:
    tuning = np.logspace(-2, 1, 10)
else:
    tuning = [1.]

implement burn in stuff here
'''

# Start loop
# parameter index
pdx = 0
# total acceptance ratio
total_accept = 0
# effective step
eff_step = 1

for i in range(0, args.iter+1):
    print ('step', i)

    # Randomly pick new parameter value if not the initial step
    if i > 0:
        current_value = parameter_space[i-1][pdx]
        # Copy previous step's info
        parameter_space[i] = parameter_space[i-1]
        # Dictionary keyword for parameter to update (only do one per iteration)
        keyword = [x for x,y in parameter_dict.items() if y==pdx]
        # Update parameter_space[i] with random proposed value
        parameter_space[i][pdx] = proposal_distribution(keyword[0],
            current_value)

    # If parameter is out of prior bounds, do not calculate likelihood.
    # Set posterior equal to the proposal (i.e., -inf).
    # If parameter is within prior bounds, calculate posterior from peak flux.
    # Peak flux posterior will return -inf if the none of the random GRBs are
    # "detected."
    if np.isinf(parameter_space[i][pdx]) == True:
        parameter_space[i][post_idx] = parameter_space[i][pdx]
    else:
        parameter_space[i][post_idx] = Simulation(parameter_space[i],
            parameter_dict, num_param, detector=detector,
            redshifts=redshifts, luminosities=luminosities, obs_pf=obs_pf,
            dl=dl, options=options, vol_arr=v_comov, kc=kcorr, plot_GRB=plot,
            prior=args.prior)

    if i > 0:
        # Previous parameter_space entry with parameter pdx
        step_back = np.maximum(0,i-num_param)

        # Metro-Hastings decision algorithm
        # Change only if parameter is rejected, which should be most of the time
        # (probably more efficient if only change if parameter is accepted,
        # but doesn't work with my code right now)

        # Give proposed parameter, LL with proposed parameter, previous
        # parameter acceptance ratio, and total acceptance ratio
        parameter_space[i][pdx], parameter_space[i][post_idx], \
        parameter_space[i][p_accept_idx], total_accept = \
            metro_hastings(current_param=parameter_space[i-1][pdx],
                proposed_param=parameter_space[i][pdx],
                current_post=parameter_space[i-1][post_idx],
                proposed_post=parameter_space[i][post_idx],
                total_accept=total_accept,
                p_accept=parameter_space[step_back][p_accept_idx])

        # Record the effective step of the parameter
        parameter_space[i][eff_step_idx] = eff_step

    if i > 0:
        # Restart parameter loop if reached the end
        pdx += 1
        if pdx > num_param-1:
            pdx = 0
            eff_step += 1

    # Every so often, save the results
    if np.mod(i, 10000) == 0:
        if args.filename is not None:
            print ('Saving '+file+'...')
            np.save(file, parameter_space)
print (parameter_space[i])

print("--- Runtime of %s seconds ---\n" % (time.time() - sim_time))
print ('Acception fraction:', total_accept/(args.iter))


# Save file of parameters, LLRs, and diagnostics
if args.filename is not None:
    np.save(file, parameter_space)
    print (file)
