import glob
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from math_functions import autocorrelation

def plot_samples(x, pdfx, pdfy, xlog=False, ylog=False, xlabel=None, ylabel=None):
    if xlog is not False:
        pdfy = np.array(pdfy)
        pdfy = pdfy / pdfy[0]
        plt.plot(pdfx, pdfy, label='pdf', color='red')

        plt.xscale('log')
        bins = np.logspace(50, 54, 30)
        yvals, histbins = np.histogram(x, bins=bins)
        bindiff = np.diff(histbins)
        print (np.diff(histbins))
        print (np.diff(np.log10(histbins)))
        newy = yvals / yvals[0]
        #newy = yvals / (len(x) * np.diff(histbins))
        plt.bar(histbins[:-1], newy, width=bindiff, align='edge', edgecolor='C0', color='None', label='samples')
    elif xlog is False:
        plt.hist(x, bins=30, histtype='step', density=True)
    if ylog is not False:
        plt.yscale('log')
    else:
        plt.plot(pdfx, pdfy, color='red', label='pdf')
    plt.legend()
    plt.show()
    plt.close()
    return

def plot_rate(redshifts, rate):
    fig = plt.figure(constrained_layout=True)
    grid = plt.GridSpec(1, 4, figure=fig)
    first = plt.subplot(grid[0, :2])
    first.plot(redshifts, rate)
    plt.xlabel('Redshift')
    plt.ylabel(r'$R_{GRB;obs}$ [$Gpc^{-1} yr^{-1}$]')
    plt.title('Obs GRB Rate')
    '''
    # Intrinsic comoving grb rate [$Gpc^{-3} yr^{-1}$]
    sum = quad(int_grb_rate, redshifts[0], redshifts[-1], args=(z_star, n1, n2, rho0))
    pdf = [int_grb_rate(z, r_0=rho0/sum[0], z_star=z_star, n1=n1, n2=n2) for z in redshifts]
    second = plt.subplot(grid[0, 2:])
    second.plot(redshifts, pdf)
    plt.xlabel('Redshift')
    plt.ylabel(r'$R_{GRB}$ [$Gpc^{-1} yr^{-1}$]')
    plt.title('Normalized Comoving GRB Rate')
    '''
    plt.show()
    plt.close()
    return

def plot_cdf(xs, cdf_normed, xlabel, xlog=False, ylog=False):
   with plt.xkcd():
        if xlog is not False:
            plt.xscale('log')
        if ylog is not False:
            plt.yscale('log')
        plt.plot(xs, cdf_normed)
        plt.xlabel(xlabel)
        plt.ylabel('CDF')
        plt.grid()
        plt.show()
        plt.close()
        return

def plot_det_thresh(s=1, w=1):
    # Tie logistic function to peak flux distribution
    y = np.logspace(-1, 1, 100)
    b = (s*y**w) / (1 + s*y**w)
    plt.xscale('log')
    plt.plot(y,b)
    plt.grid()
    #props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    #plt.text(0.2,0.9, r'y = $\frac{Ax^B}{(1+Ax^B)}$', fontsize=14)
    plt.show()
    return

def plot_peak_flux_color(L, z, corr_pf):
    plt.scatter(L, z, c=corr_pf, cmap='hsv', norm=matplotlib.colors.LogNorm())
    plt.xscale('log')
    #plt.yscale('log')
    plt.xlabel('Luminosity')
    plt.ylabel('Redshift')
    plt.colorbar(label='Peak flux')
    plt.tight_layout()
    plt.show()
    plt.close()
    return

def plot_peak_flux_norm(model_counts, data_counts, bins):
    model_norm = np.array([n/np.sum(model_counts) for n in model_counts])
    data_norm = np.array([n/np.sum(data_counts) for n in data_counts])
    plt.step(bins[:-1], model_norm, where='post', label='Model (All)')
    plt.step(bins[:-1], data_norm, where='post', label='Data (All)')
    plt.xscale('log')
    plt.legend()
    plt.show()
    plt.close()
    return

def plot_peak_flux(coll_all, merg_all, coll_mod, merg_mod, data,
            coll_model_label, merg_model_label, total_cut_model):
    total_model = np.concatenate([coll_all,merg_all])
    bins=np.logspace(-1, 3, 70)
    #bins = np.logspace(-10, -5, 200)
    # Unnormalized plots of all data
    plt.hist(coll_all, bins=bins, histtype='step', label=coll_model_label)
    plt.hist(merg_all, bins=bins, histtype='step', label=merg_model_label)
    plt.hist(total_model, bins=bins, histtype='step', label='All Model')
    plt.hist(data, bins=bins, histtype='step', linewidth=2, label='Data')
    plt.hist(total_cut_model, bins=bins, histtype='step', linewidth=2, label='Model Total Detected')
    plt.xscale('log')
    plt.xlabel('1s Peak Flux [1-10,000 keV]')
    plt.ylabel('Number of GRBs')
    plt.legend()
    #plt.show()
    plt.close()

    # Unnormalized cut peak flux
    bins=np.logspace(-1, 3, 40)
    #bins = np.logspace(-11, -7, 200)
    #plt.hist(coll_mod, bins=bins, histtype='step', label=coll_model_label)
    #plt.hist(merg_mod, bins=bins, histtype='step', label=merg_model_label)
    plt.hist(total_cut_model, bins=bins, histtype='step', label='Model (All)')
    plt.hist(data, bins=bins, histtype='step', label='Data (All)')
    plt.xscale('log')
    plt.legend()
    plt.show()
    plt.close()
    return


# Obs luminosities
def plot_obs_lum(final_l):
    bins = np.logspace(50,54,100)
    plt.hist(final_l, bins=bins)
    plt.xscale('log')
    plt.show()
    plt.close()
    return


# obs grb rate
def plot_obs_rate(final_z):
    plt.hist(final_z, bins=30)
    plt.show()
    plt.close()
    return


def trace(step, x, label=None):
    # Trace plots of parameters
    plt.plot(step, x, label=label)
    plt.legend()
    plt.xlabel('iteration')
    plt.title('Trace Plot')
    plt.show()
    plt.close()
    return


def plot_autocorrelation(x, label=None):
    # Compute autocoreelation function for chain
    lags = np.arange(0, 1000)
    rho = autocorrelation(x, lags)
    plt.plot(lags, rho, label=label)
    plt.xlabel(r'lag $k$')
    plt.ylabel(r'estimated autocorrelation $\rho(k)$')
    plt.title('Correlation between every kth sample of MCMC chain')
    plt.legend()
    plt.show()
    plt.close()


def llr(i, x):
    # Plot LLR as a function of step
    plt.plot(i, x)
    plt.xlabel('iteration')
    plt.ylabel('posterior value')
    plt.yscale('log')
    plt.show()
    plt.close()
    return


def plot_parameter_grid(first, second=None, third=None, fourth=None):
    fig = plt.figure(figsize=(9,9))
    grid = plt.GridSpec(4, 4)

    # LEFT-MOST ROW
    # Histogram of first parameter
    f = plt.subplot(grid[0, 0])
    h = f.hist(first, bins=30, histtype='step')
    plt.vlines(np.median(first), 0, np.max(h[0]), linestyle='--', color='red')
    plt.ylabel('coll z1')
    # Scatter plot of first parameter and second parameter
    f_corr = plt.subplot(grid[1, 0])
    f_corr.scatter(first, second, s=2)
    plt.ylabel('coll z2')
    # Scatter plot of first parameter and third parameter
    f_sec_corr = plt.subplot(grid[2, 0])
    f_sec_corr.scatter(first, third, s=2)
    plt.ylabel('coll z*')
    # Scatter plot of first parameter and fourth parameter
    f_thi_corr = plt.subplot(grid[3, 0])
    f_thi_corr.scatter(first, fourth, s=2)
    plt.ylabel('coll rho0')
    plt.xlabel('coll z1')

    # MIDDLE ROW
    # Histogram of second parameter
    s = plt.subplot(grid[1, 1])
    h = s.hist(second, bins=20, histtype='step')
    plt.vlines(np.median(second), 0, np.max(h[0]), linestyle='--', color='red')
    # Scatter plot of second parameter and third
    second_corr = plt.subplot(grid[2, 1])
    second_corr.scatter(second, third, s=2)
    #plt.yticks([])
    plt.xlabel('coll z2')

    # Third row
    t = plt.subplot(grid[2, 2])
    h = t.hist(third, bins=20, histtype='step')
    plt.vlines(np.median(third), 0, np.max(h[0]), linestyle='--', color='red')
    #plt.yticks([])
    plt.xlabel('coll z*')

    # Finish
    plt.tight_layout()
    plt.show()
    plt.close()



# Read file
def plot_results(filename):
    results = np.load(filename)
    parameters = {
        0: "coll_z1",
        1: "coll_z2",
        2: "coll_z*",
    }
    i = np.arange(1, len(results))
    z1 = results[:,0][1:]
    z2 = results[:,1][1:]
    zstart = results[:,2][1:]
    LLR = results[:,6][1:]
    accepted = results[:,7][1:]
    effective_step = results[:,8][1:]
    acceptance_ratio = accepted / effective_step
    '''
    # Plot "running mean plots"
    # calculate the mean as a function of iteration
    print ('check running mean plot')
    summed =  np.cumsum(z1)
    run_mean = [summed[i]/float(i+1) for i in range(len(summed))]
    plt.plot(i, run_mean, label='z1')
    summed =  np.cumsum(z2)
    run_mean = [summed[i]/float(i+1) for i in range(len(summed))]
    plt.plot(i, run_mean, label='z2')
    summed =  np.cumsum(zstart)
    run_mean = [summed[i]/float(i+1) for i in range(len(summed))]
    plt.plot(i, run_mean, label='z*')
    plt.legend()
    plt.xlabel('iteration')
    plt.title('Running Mean Plot')
    plt.show()
    plt.close()
    '''
    return

def track_acceptance(file):
    results = np.load(file)
    accepted = results[:,7][1:]
    effective_step = results[:,8][1:]
    return accepted / effective_step








'''
tuning = tuning = np.logspace(-2, 1, 10)
files = ['results/test_tuning_0.01.npy', 'results/test_tuning_0.02.npy',
        'results/test_tuning_0.04.npy', 'results/test_tuning_0.1.npy',
        'results/test_tuning_0.2.npy', 'results/test_tuning_0.4.npy',
        'results/test_tuning_1.0.npy', 'results/test_tuning_2.npy',
        'results/test_tuning_4.npy', 'results/test_tuning_10.0.npy']

z1s = []
z2s = []
zstars = []
for f in range(len(files)):
    ratios = track_acceptance(files[f])
    z1s.append(np.median(ratios[::3]))
    z2s.append(np.median(ratios[1::3]))
    zstars.append(np.median(ratios[2::3]))

plt.plot(tuning, z1s, label='z1')
plt.scatter(tuning, z1s)
plt.plot(tuning, z2s, label='z2')
plt.scatter(tuning, z2s)
plt.plot(tuning, zstars, label='z2')
plt.scatter(tuning, zstars)
plt.xscale('log')
plt.legend()
#plt.show()
plt.close()
'''

'''
if options.getboolean('testing') is not False:
    plt.plot(redshifts, coll_redshift_pdf, color='red')
    bins = np.linspace(0,10,40)
    n = np.histogram(coll_redshift_sample, bins=bins)
    normalization = np.sum(n[0] * np.diff(n[1]))
    weights = [1./normalization] * len(coll_redshift_sample)
    c = plt.hist(coll_redshift_sample, bins=bins, weights=weights)
    plt.xlabel('Redshift')
    plt.ylabel('Probability Density')
    plt.grid()
    plt.show()
    plt.close()

    # Cumulative distribution of samples
    plt.plot(redshifts, coll_redshift_cdf, color='red')
    weights = [1/np.sum(n[0])] * len(coll_redshift_sample)
    plt.hist(coll_redshift_sample, bins=bins, cumulative=True, weights=weights)
    plt.grid()
    plt.xlabel('Redshift')
    plt.ylabel('cdf')
    plt.show()
    plt.close()
'''
