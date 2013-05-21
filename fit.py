#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
import logging
import correlator
import build_corr
import pylab

from scipy import linalg
from scipy.optimize import leastsq
from scipy.optimize import fmin
# from scipy.optimize import minimize
Nt = 128

def single_exp(v, x):
    return (v[1] * np.exp((-1.0) * v[0] * x))

def periodic_exp(v, x):
    return (v[1] * (np.exp((-1.0) * v[0] * x) + np.exp(v[0] * (x-(Nt/2.0)) )))


def cosh(v, x):
    #return ((2*v[1])/np.exp(v[0]*Nt/2.0) * np.cosh((-1.0)* v[0]*((x-(Nt/2.0)))))
    return (v[1] * np.cosh((-1.0)* v[0]*((x-(Nt/2.0)))))

def fit(fn, cor, tmin, tmax):
    fun = lambda v, mx, my: (fn(v, mx) - my)

    initial_guess = [0.01, 1.0]
    x = np.arange(32)[tmin:tmax]
    ave_cor = cor.average_sub_vev()
    y = [ave_cor[t] for t in range(tmin, tmax)]
    original_ensamble_params, success = leastsq(fun, initial_guess, args=(x, y), maxfev=100000)
    if not success:
        raise ValueError()
    initial_guess = original_ensamble_params

    def cov_fit(correlator,guess):
        ave_cor = i.average_sub_vev()
        y = [ave_cor[t] for t in range(tmin, tmax)]
        cov = covariance_matrix(i, tmin, tmax)

        inv_cov = bestInverse(cov)

        cov_fun = lambda v: np.sum(((ave_cor[t] - fn(v, t)) * inv_cov[t - tmin][tp - tmin] * (ave_cor[tp] - fn(v, tp))
                                    for t in range(tmin, tmax) for tp in range(tmin, tmax)))
        uncorrelated_fit_values, success = leastsq(fun, guess, args=(x, y), maxfev=100000)
        if not success:
            raise ValueError()

        guess = uncorrelated_fit_values
        #would like to use minimize, but seems to be not installed
        results = fmin(cov_fun, guess, ftol=1.E-7, maxfun=1000000, full_output=1, disp=0, retall=0)
        covariant_fit, fit_info, flag = results[0], results[1:-1], results[-1]
        # covariant_fit = minimize(cov_fun, initial_guess)
        logging.debug("Fit results: f() ={}, Iterations={}, Function evaluations={}".format(*fit_info))
        if flag:
            logging.error("Fitter flag set to {}. Error!")
            raise RuntimeError("Fitter failed")

        logging.debug("Covariant fit {}, regular fit {}".format(repr(covariant_fit),
                                                                repr(uncorrelated_fit_values)))

        return(covariant_fit)

    boot_masses = []
    boot_amps = []
    for i in bootstrap_ensamble(cor):
        fitted_params = cov_fit(i, original_ensamble_params)
        boot_masses.append(fitted_params[0])
        boot_amps.append(fitted_params[1])

    original_ensamble_correlatedfit = cov_fit(cor, original_ensamble_params)

    print 'Uncorelated total fit: ', original_ensamble_params
    print 'Correlated total fit:  ', original_ensamble_correlatedfit
    print 'boot parameter masses: ', np.mean(boot_masses), np.var(boot_masses)

    v = (np.mean(boot_masses), np.mean(boot_amps))
    cov = covariance_matrix(cor, tmin, tmax)
    inv_cov = bestInverse(cov)
    chi_sqr = np.sum(((ave_cor[t] - fn(v, t)) * inv_cov[t - tmin][tp - tmin] * (ave_cor[tp] - fn(v, tp))
                      for t in range(tmin, tmax) for tp in range(tmin, tmax)))

    print "Chi_Sq", chi_sqr
    dof = len(x) - 2
    print "Chi_Sq / dof", chi_sqr/dof

    return (v, boot_masses)


def plot_fit(fn, cor, tmin, tmax, filename=None):
    X = np.linspace(tmin, tmax, 200 * 5)
    fitted_params, boot_masses = fit(fn, cor, tmin, tmax)

    plt.figure()
    corplot = plt.subplot(211)

    corplot.errorbar(cor.times, cor.average_sub_vev().values(), yerr=cor.jackknifed_errors(), fmt='o')
    corplot.plot(cor.times, cor.average_sub_vev().values(), 'ro', X, fn(fitted_params, X))
    plt.ylim([0, max(cor.average_sub_vev().values())])
    emass = cor.effective_mass(1)
    emass_errors = cor.effective_mass_errors(1).values()
    emassplot = plt.subplot(212)
    dataplt = emassplot.errorbar(emass.keys(), emass.values(), yerr=emass_errors, fmt='o')
    fitplt = emassplot.errorbar(X, np.mean(boot_masses) * np.ones_like(X), yerr=np.std(boot_masses))
    plt.legend([dataplt, fitplt], ["data", u"fit m={:.5f}\xb1{:.5f}".format(np.mean(boot_masses), np.std(boot_masses))])
    plt.ylim([0, 1])
    plt.xlim([0, tmax + 2])
    if(filename):
        plt.savefig(filename)
    else:
        plt.show()


def make_fake_cor():
    cfgs = list(range(50))
    times = list(range(16))
    data = {}
    vev = {}
    for c in cfgs:
        vev[c] = 0.0
        tmp = {}
        for t in times:
            tmp[t] = 5.0 * np.exp((-0.5) * t) + 6.0 * np.exp((-1.5) * t)
            tmp[t] += (pylab.rand() * 0.001) * tmp[t]
        data[c] = tmp

    return correlator.Correlator.fromDataDicts(data, vev, vev)


def bootstrap_cfgs(cor):
    return np.random.choice(cor.configs, size=len(cor.configs))


def bootstrap(cor):
    newcfgs = bootstrap_cfgs(cor)
    newcor = {i: cor.get(config=c) for i, c in enumerate(newcfgs)}
    # print cor.vev1
    newvev1 = {i: cor.vev1[c] for i, c in enumerate(newcfgs)}
    newvev2 = {i: cor.vev2[c] for i, c in enumerate(newcfgs)}
    return correlator.Correlator.fromDataDicts(newcor, newvev1, newvev2)


def bootstrap_ensamble(cor, N=256):
    return [bootstrap(cor) for i in range(N)]


def covariance_matrix(cor, tmin, tmax):
    # print(cor.data)
    # reformat = np.array([c.values() for c in cor.data.values()])
    # print(reformat)
    nm1 = (1.0 / (len(cor.configs) - 1))
    mymat = np.identity(tmax - tmin)
    aoc = cor.average_over_configs()
    for i in range(tmin, tmax):
        ave_i = aoc[i]
        for j in range(i, tmax):
            ave_j = aoc[j]
            mymat[i - tmin][j - tmin] = np.mean([(cor[n][i] - ave_i) * (cor[n][j] - ave_j)
                                                 for n in cor.configs]) * nm1
            mymat[j - tmin][i - tmin] = mymat[i - tmin][j - tmin]

    return mymat


def CholeskyInverse(t):
    """
    Computes inverse of matrix given its Cholesky upper Triangular
    decomposition t.
    """
    nrows = len(t)
    B = np.zeros_like(t)

    # Backward step for inverse.
    for j in reversed(range(nrows)):
        tjj = t[j][j]
        S = np.sum([t[j][k] * B[j][k] for k in range(j + 1, nrows)])
        B[j][j] = 1.0 / tjj**2 - S / tjj  # flake8: noqa
        for i in reversed(range(j)):
            B[j][i] = B[i][j] = -np.sum([t[i][k] * B[k][j] for k in range(i + 1, nrows)]) / t[i][i]
    return B

class InversionError(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


def bestInverse(M):
    TOLERANCE = 1.E-5

    def invert_error(i):
        return np.max(np.abs((np.dot(M, i) - np.identity(len(i)))))

    inv = linalg.inv(M)         # Will raise error if singular
    error = invert_error(inv)

    try:
        chol = linalg.cholesky(M, check_finite=False )
    except np.linalg.linalg.LinAlgError:
        logging.error("Not positive definite!")

    else:
        chol_inv = CholeskyInverse(chol)
        chol_error = invert_error(chol_inv)

        logging.debug("Inversions errors were inv={}, chol={}".format(error, chol_error))

        if chol_error > TOLERANCE and error > TOLERANCE:
            raise InversionError("Could not invert within tolerance")

        if chol_error < error:
            logging.debug("Using choleskey inverse")
            inv = chol_inv
        else:
            logging.debug("Using standard inverse")
    return inv


if __name__ == "__main__":
    logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.DEBUG)

    corrfile = "/home/bfahy/r2/testdata/myformatcorsnk-etap000DDL7Egp1_src-etap000DDL7Egp1.dat"
    # srcvevfile = snkvevfile = "/home/bfahy/r2/combined_results/lattice_26_beta_6.0/total_both/binned_500_A1++_1.looptype3_opnum0_size4_A1++_1.looptype3_opnum0_size4.vev1"

    cor = build_corr.corr_and_vev_from_files(corrfile, None, None, cfgs=50, ts=18)

    plot_fit(single_exp, cor, 10, 18)
