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
"coll mu": 12,
"coll sigma": 13,
"merg mu": 14,
"merg sigma": 15,
"posterior": 16,
"acceptance num": 17,
"effective step": 18,
}


def plot_trace(results, steps, num_param):
    for par in range(num_param):
        keyword = [x for x,y in dict.items() if y==par]
        trace(steps, results[:, par], label=keyword[0], ylog=True)
    return


num_param = 4
# Load results
filename = '../results/test_prior_log4.npy'
results = np.load(filename)

'''
#results = np.delete(results, slice(50001,100001), 0)
#results = results[1:] # cut off initial row
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

# Clean parameters
collz1 = results[:,dict["coll z1"]][::num_param]
collz2 = results[:,dict["coll z2"]][::num_param]
collzpeak = results[:,dict["coll z*"]][::num_param]
coll_local = results[:,dict["coll rho0"]][::num_param]

mergz1 = results[:,dict["merg z1"]][::num_param]
mergz2 = results[:,dict["merg z2"]][::num_param]
mergzpeak = results[:,dict["merg z*"]][::num_param]
merg_local = results[:,dict["merg rho0"]][::num_param]

coll_alpha = results[:,dict["coll alpha"]][::num_param]
coll_beta = results[:,dict["coll beta"]][::num_param]
merg_alpha = results[:,dict["merg alpha"]][::num_param]
merg_beta = results[:,dict["merg beta"]][::num_param]

#LLR = results[:,dict["posterior"]]
#accepted = results[:,dict["acceptance ratio"]]
#effective_step = results[:,dict["effective step"]]


# Show trace plots
plot_trace(results, steps, num_param)

'''
# Plot autocorrelation
plot_autocorrelation([collz1, collz2, collzpeak, coll_local], maxlag=300, label=['collapsar low-z index', 'collapsar high-z index', 'collapsar peak z', 'collapsar local rate'])
plot_autocorrelation([mergz1, mergz2, mergzpeak, merg_local], maxlag=300, label=['merger low-z index', 'merger high-z index', 'merger peak z', 'merger local rate'])

# Find autocorrelation time
print ('Input max lag:')
atime = int(input())

# Plot parameter grid (using knowledge from autocorrelation time)
labels = ['coll z1', 'coll z2', 'coll z*', 'coll rho0', 'merg z1', 'merg z2', 'merg z*', 'merg rho0']
parameter_vector = [collz1[::atime], collz2[::atime], collzpeak[::atime], coll_local[::atime], mergz1[::atime], mergz2[::atime], mergzpeak[::atime], merg_local[::atime]]
plot_grid_general(parameter_vector, labels)
'''

# Calculate correlation ?


'''
# find acceptance ratio for each parameter
#acceptance_ratio = accepted / effective_step

# Plot posterior as a function of step (should I start plotting LLR?)
# including autocorrelation?
#steps = np.arange(0, len(results))
#llr(steps, LLR)
'''
