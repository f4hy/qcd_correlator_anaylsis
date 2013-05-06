#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
import logging
import correlator
import build_corr
import pylab

from scipy.optimize import leastsq


def single_exp(v, x):
    return (v[1] * np.exp((-1.0) * v[0] * x))


def fit(fn, cor, tmin, tmax):
    # fun = lambda v, x, y: ((v[1] * np.exp((-1.0)*v[0]*x)) -y)
    # fn = lambda v, x: ((v[1] * np.exp((-1.0)*v[0]*x)) )

    fun = lambda v, x, y: (fn(v, x) - y)

    initial_guess = [2.0, 2.0]
    x = np.arange(32)[tmin:tmax]
    y = np.exp((-1.0) * np.arange(32))
    ave_cor = cor.average_sub_vev()
    y = [ave_cor[t] for t in range(tmin, tmax)]

    v, success = leastsq(fun, initial_guess, args=(x, y), maxfev=10000)
    if not success:
        raise ValueError()
    boot_masses = []
    for i in bootstrap_ensamble(cor, 64):
        ave_cor = i.average_sub_vev()
        y = [ave_cor[t] for t in range(tmin, tmax)]
        fit, success = leastsq(fun, initial_guess, args=(x, y), maxfev=10000)
        if not success:
            raise ValueError()
        boot_masses.append(fit[0])

    if not success:
        raise ValueError()
    # print "boot_masses", boot_masses
    print 'Estimater parameters: ', v
    print 'boot parameter masses: ', np.mean(boot_masses), np.var(boot_masses)

    return (v, boot_masses)


def plot_fit(fn, cor, tmin, tmax, filename=None):
    # fun = lambda v, x, y: (fn(v,x) -y)
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
    times = list(range(32))
    data = {}
    vev = {}
    for c in cfgs:
        vev[c] = 0.0
        tmp = {}
        for t in times:

            tmp[t] = 5.0 * np.exp((-0.5) * t) + 6.0 * np.exp((-1.5) * t)
            tmp[t] *= 1 + (pylab.rand() * 0.001)
        # print tmp
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


def bootstrap_ensamble(cor, N=64):
    return [bootstrap(cor) for i in range(N)]


def covariance_matrix(cor):
    pass

if __name__ == "__main__":
    # pylab.plot(np.arange(32), make_fake_cor()[0])
    # pylab.show()
    # cor = make_fake_cor()
    # ave_cor = cor.average_sub_vev()
    # print bootstrap(cor)

    # print ave_cor
    # print bootstrap(cor).average_over_configs()

    # print bootstrap_ensamble(cor)
    corrfile = "/home/bfahy/r2/combined_results/lattice_26_beta_6.0/total_both/binned_500_A1++_1.looptype3_opnum0_size4_A1++_1.looptype3_opnum0_size4.cor"
    srcvevfile = snkvevfile = "/home/bfahy/r2/combined_results/lattice_26_beta_6.0/total_both/binned_500_A1++_1.looptype3_opnum0_size4_A1++_1.looptype3_opnum0_size4.vev1"
    cor = build_corr.corr_and_vev_from_files(corrfile, srcvevfile, snkvevfile)

    plot_fit(single_exp, cor, 2, 20)
