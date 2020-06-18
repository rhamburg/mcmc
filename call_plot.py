from plot import *

dict = {
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
    "acceptance ratio": 13,
    "effective step": 14,
    }

num_param = 12
# Load results
filename = 'results/test_luminosity_prior.npy'
results = np.load(filename)
#results = np.delete(results, slice(50001,100001), 0)
#results = results[1:] # cut off initial row

'''
# Optional concatenating of second array
file2 = '/Users/rachelhamburg/Dropbox/FermiWork Team Folder/GradSchool/Dissertation/project/mcmc/results/test8b.npy'
results2 = np.load(file2)
results2 = np.delete(results2, slice(50001,100001), 0)
results2 = results2[1:] # cut off the initial row, which was last of previous
results = np.concatenate([results, results2])
'''

# Parameters
steps = np.arange(0, len(results))
tot_eff_steps = int(len(results) / num_param)
i = np.arange(0, tot_eff_steps)

z1 = results[:,dict["coll z1"]][::num_param]
z2 = results[:,dict["coll z2"]][::num_param]
zpeak = results[:,dict["coll z*"]][::num_param]
coll_local = results[:,dict["coll rho0"]][::num_param]

LLR = results[:,dict["posterior"]]
accepted = results[:,dict["acceptance ratio"]]
effective_step = results[:,dict["effective step"]]

# Plot trace plot

for par in range(num_param):
    keyword = [x for x,y in dict.items() if y==par]
    trace(steps, results[:, par], label=keyword[0])
'''
trace(i, z1, label='collapsar z1')
trace(i, z2, label='collapsar z2')
trace(i, zpeak, label='collapsar z*')
trace(i, coll_local, label='collapsar rho0')
'''

# Plot autocorrelation
'''
plot_autocorrelation([z1, z2, zpeak, local], maxlag=300,
    label=['collapsar low-z index', 'collapsar high-z index',
    'collapsar peak z', 'collapsar local rate'])
'''


# Find autocorrelation time
# atime

'''
# Plot parameter grid (using knowledge from autocorrelation time)
atime = 100
z1 = z1[::atime]
z2 = z2[::atime]
zpeak = zpeak[::atime]
local = local[::atime]
print (len(z1))
plot_parameter_grid(z1, second=z2, third=zpeak, fourth=local)
'''



# find acceptance ratio for each parameter
#acceptance_ratio = accepted / effective_step

# Plot posterior as a function of step (should I start plotting LLR?)
# including autocorrelation?
#steps = np.arange(0, len(results))
#llr(steps, LLR)
