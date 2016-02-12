#!/usr/bin/env python
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import logging
import correlator
import build_corr
import argparse
import os
import math

from parser_fit import fitparser, functions
from fit_parents import InvalidFit
from copy import deepcopy

from scipy import linalg
from scipy import stats
from scipy.special import gammaincc
from scipy.optimize import leastsq
# from scipy.optimize import fmin
# from scipy.optimize import fmin_slsqp
# from scipy.optimize import fmin_l_bfgs_b
# from scipy.optimize import minimize

import progress_bar

OUTPUT = 25
ALWAYSINFO = 26
logging.addLevelName(OUTPUT, "OUTPUT")
logging.addLevelName(ALWAYSINFO, "INFO")

Nt = 128
NBOOTSTRAPS = 1000

def fit(fn, cor, tmin, tmax, filestub=None, bootstraps=NBOOTSTRAPS, return_quality=False,
        return_chi=False, writecor=True, options=None):
    if(tmax-tmin < len(fn.parameter_names)):
        raise InvalidFit("Can not fit to less points than parameters")

    if options.random:
        logging.info("Setting random seed to %s", options.random)
        np.random.seed(options.random)

    results = logging.getLogger("results")
    if filestub and not results.handlers:
        filename = filestub+".stats"
        filehandler = logging.FileHandler(filename)
        filehandler.level = OUTPUT
        results.addHandler(filehandler)
        logging.info("Writing output to file {}".format(filename))

    results.info("Fitting data to {} from t={} to t={} using {} bootstrap samples".format(
        fn.description, tmin, tmax, bootstraps))

    #tmax = tmax+1  # I use ranges, so this needs to be offset by one
    fitrange = range(tmin, tmax+1)
    fun = lambda v, mx, my: (fn.formula(v, mx) - my)
    initial_guess = fn.starting_guess(cor, options.period, tmax, tmin)
    logging.info("Starting with initial_guess: {}".format(repr(initial_guess)))

    x = np.array(fitrange)
    dof = len(x) - len(fn.parameter_names)
    orig_ave_cor = cor.average_sub_vev()
    y = [orig_ave_cor[t] for t in fitrange]

    if fn.subtract:
        logging.debug("before suctracted correlator is:")
        logging.debug(cor.average_sub_vev())
        ccor = deepcopy(cor)
        cor = ccor
        fn.subtract = tmin - 1
        cor.subtract(tmin-1)
        logging.debug("subtracted correlator is:")
        logging.debug(cor.average_sub_vev())

    orig_ave_cor = cor.average_sub_vev()
    y = [orig_ave_cor[t] for t in fitrange]
    original_ensamble_params, success = leastsq(fun, initial_guess, args=(x, y), maxfev=10000)
    if options.debugguess:
        #return original_ensamble_params, [0.01, 0.01, 0.01, 0.01] # For testing initila guess in plot
        if options.plot:
            plot_fit(fn, cor, tmin, tmax, options, initial_guess)
        return initial_guess, [0.01, 0.01, 0.01, 0.01]  # For testing initila guess in plot
    if not success:
        raise InvalidFit("original exnamble leastsq failed")
    if options.first_pass:
        initial_guess = original_ensamble_params
        logging.info("initial_guess after first pass: {}".format(repr(initial_guess)))

    def cov_fit(correlator, guess):
        ave_cor = correlator.average_sub_vev()
        y = [ave_cor[t] for t in fitrange]


        if options.debug_uncorrelated:
            logging.debug("Using uncorrlated")
            # cov = covariance_matrix(correlator, tmin, tmax+1)
            # cov = np.diag(np.diag(cov))
            jke = correlator.jackknifed_errors()
            cov = np.diag([jke[t]**2 for t in fitrange])
        elif options.debug_jkcov:
            cov = jk_covariance_matrix(correlator, tmin, tmax+1)
        else:
            cov = covariance_matrix(correlator, tmin, tmax+1)

        inv_cov = bestInverse(cov)


        if options.debug_identcov:
            results.log(30, "using identcov debug option")
            inv_cov = np.identity(len(cov))


        aoc = np.array([ave_cor[t] for t in fitrange])
        #logging.debug("guess {}".format(str(guess)))

        def cov_fun(g):
            """ Function to be minizied. computed using matrix mult"""
            vect = aoc - fn.formula(g, x)
            return vect.dot(inv_cov).dot(vect)
        if options.first_pass:
            uncorrelated_fit_values, success = leastsq(fun, guess, args=(x, y), maxfev=100000)
            if not success:
                raise InvalidFit("leastsq failed")
            logging.debug("firstpass guess {}".format(str(uncorrelated_fit_values)))
            if guess[0] < 0.0:
                logging.warn("first pass found mass to be negative {}, lets not use it".format(guess[0]))
            else:
                guess = uncorrelated_fit_values

            if len(guess) > 2 and guess[2] < 0.0:
                logging.warn("first pass found mass2 to be negative {}, lets flip it".format(guess[2]))
                logging.info("first pass results are {}".format(repr(guess)))
                guess[2] = -guess[2]

        def clamp(n, minn, maxn):
                return max(min(maxn, n), minn)
        bounded_guess = [clamp(g, b[0], b[1]) for g, b in zip(guess, fn.bounds)]
        #logging.debug("guess {}, bounded guess {}".format(repr(guess), repr(bounded_guess)))

        m = fn.custom_minuit(aoc, inv_cov, x, guess=bounded_guess)
        #m.set_strategy(2)
        migradinfo = m.migrad()
        minuit_results = [m.values[name] for name in fn.parameter_names]
        chisqr = migradinfo[0]["fval"]
        if m.get_fmin().is_valid:
            return minuit_results, chisqr
        else:
            logging.error("minuit failed!!")
            logging.error("was at {}".format(minuit_results))
            raise InvalidFit("minuit failed")

    # end cov_fit

    original_ensamble_correlatedfit, original_ensamble_chisqr = cov_fit(cor, initial_guess)
    isvalidfit = fn.valid(original_ensamble_correlatedfit)
    if not isvalidfit:
        raise InvalidFit("Full ensamble failed")


    boot_params = []
    boot_chisqr = []
    failcount = 0
    attempted = 0

    for strap in bootstrap_ensamble(cor, N=bootstraps, filelog=filestub, jackknife=options.jackknife):

        if options.jackknife:
            bootstraps = len(cor.configs)

        pb = progress_bar.progress_bar(bootstraps)

        attempted +=1
        pb.update(attempted)
        if options.reguess:
            newguess = fn.starting_guess(strap, options.period, tmax, tmin)
        else:
            newguess = initial_guess
        fitted_params, fitted_chisqr = cov_fit(strap, newguess)
        if fitted_params is not None:
            boot_params.append(fitted_params)
            boot_chisqr.append(fitted_chisqr)
            logging.debug("bootstrap converged")
            if options.write_each_boot:
                write_fitted_cor(fn, strap, tmin, tmax, options, fitted_params, postfix=".bootstrap{}".format(attempted))
            if options.debug:
                plot_fit(fn, strap, tmin, tmax, options, fitted_params, postfix=".bootstrap{}".format(attempted))
        else:
            logging.error("bootstrap failed to converge!")
            raise InvalidFit("one bootstrap failed")
            #raw_input("test")
            failcount+=1
            logging.debug("fails:{} attempts:{}, ratio:{}".format(failcount, attempted, failcount/float(attempted)))
            if failcount/float(attempted) > 0.15 and attempted > 40:
                raise InvalidFit("more than 20% of boostraps failed to converge")
        del strap
    pb.done()


    if failcount > 0:
        logging.warn("{} bootstraps did not converge!".format(bootstraps-len(boot_params)))
    if len(boot_params) < bootstraps * 0.9:
        logging.error("More that 10% of the straps failed")
        raise InvalidFit("more than 10% of boostraps failed to converge")

    if options.histo:
        plot_histograms(fn.parameter_names, boot_params, options)

    results.info('')
    results.info('Uncorelated total fit: %s', {n: p for n, p in zip(fn.parameter_names, original_ensamble_params)})
    results.info('Correlated total fit:  %s', {n: p for n, p in zip(fn.parameter_names, original_ensamble_correlatedfit)})

    factor = 1
    if options.jackknife:
        factor = np.sqrt(len(cor.configs)-1)


    boot_averages = np.mean(boot_params, 0)
    boot_std = factor * np.std(boot_params, 0)
    boota = np.array(boot_params)
    upper_quartiles = [stats.scoreatpercentile(boota[:, i], 75) for i in range(len(boot_averages))]
    medians = [stats.scoreatpercentile(boota[:, i], 50) for i in range(len(boot_averages))]
    lower_quartiles = [stats.scoreatpercentile(boota[:, i], 25) for i in range(len(boot_averages))]
    inter_range = [stats.scoreatpercentile(boota[:, i], 75) - stats.scoreatpercentile(boota[:, i], 25) for i in range(len(boot_averages))]

    for name, boot, original, err in zip(fn.parameter_names, boot_averages,
                                         original_ensamble_correlatedfit, boot_std):
        bias = abs(boot-original)
        percent_bias = abs(boot-original)/original
        results.info('Bootstrap Bias in {:<10}: {:.3%}'.format(name, percent_bias))
        if bias > err*2 and (err > 0.0):
            results.error('Bootstrap Bias in {:<10}: {:.3%}'.format(name, percent_bias))
            results.error("Bootstrap average does not agree with ensamble average!"
                          "\nNot enough statistics for this for to be valid!!!\n")
            if not options.unsafe:
                results.critical("Exiting! Run with --unsafe to fit anyway")
                raise InvalidFit("Bootstrap average does not agree with ensamble average")

    for name, ave, med, std, iqr in zip(fn.parameter_names, boot_averages, medians, boot_std,
                                        inter_range):
        skew = abs(ave-med)/ave
        dist_skew = abs(std-iqr)/iqr
        if skew > 1.0:
            results.error("{}: diff of bstrap average and bstrap med is {:.3%}".format(name, skew))
            results.error("Bootstrap distrubtion is skewed!!")
            if not options.unsafe:
                results.critical("Exiting! Run with --unsafe to fit anyway")
                raise InvalidFit("Bootstrap average does not agree with ensamble average")
        else:
            results.info("{}: diff of bstrap average and bstrap med is {:.3%}".format(name, skew))

    results.info("")
    try:
        results.log(OUTPUT, "Fit ranges ({}), ({})".format(*fn.ranges))
    except Exception as e:
        pass

    results.log(OUTPUT, "Full esamble fitted parameters t=%2d to %2d---------------------", tmin, tmax)
    results.log(OUTPUT, "Name      : Average,        STD,           (1st Quart, Median, 3rd Quart, IQR)")
    for name, ave, std, low, med, up, iqr in zip(fn.parameter_names, boot_averages, boot_std,
                                                 upper_quartiles, medians, lower_quartiles, inter_range):
        results.log(OUTPUT, u"{:<10}: {:<15.10f} \u00b1 {:<10g}   ({:<9.6f}, {:<9.6f}, {:<9.6f}, {:<9.6f})".format(name, ave, std, low, med, up, iqr))
    results.log(OUTPUT, "--------------------------------------------------------")

    v = original_ensamble_correlatedfit

    chi_average = np.mean(boot_chisqr, 0)
    chi_median = np.median(boot_chisqr, 0)
    chi_min = min(boot_chisqr)
    chi_std = np.std(boot_chisqr, 0)
    chi_range = stats.scoreatpercentile(boot_chisqr, 84) - stats.scoreatpercentile(boot_chisqr, 16)
    dof = len(x) - len(fn.parameter_names)
    chi_sqr = original_ensamble_chisqr
    results.log(OUTPUT, u'\u03c7\u00b2 ={},   \u03c7\u00b2 / dof = {}, Qual {}\n'.format(
        chi_sqr, chi_sqr/dof, quality_of_fit(dof, chi_sqr)))

    logging.debug("chiave:{}, chi_med:{}, chi_min:{}, chi_std:{}, chi_range{}".format(
        chi_average, chi_median, chi_min, chi_std, chi_range))


    if bootstraps > 1 and filestub:
        bootfilename = filestub+".boot"
        if options.jackknife:
            bootfilename = filestub+".jack"
        results.info("writing each bootstrap parameter to {}".format(bootfilename))
        with open(bootfilename, 'w') as bootfile:
            str_ensamble_params = ", ".join([str(p) for p in original_ensamble_correlatedfit])
            bootfile.write("#bootstrap, {}, \t ensamble mean: {}\n".format(", ".join(fn.parameter_names), str_ensamble_params))
            for i, params in enumerate(boot_params):
                strparams = ", ".join([str(p) for p in params])
                bootfile.write("{}, {}\n".format(i, strparams))

    if options.output_stub and writecor:
        write_fitted_cor(fn, cor, tmin, tmax, options, boot_averages, errors=boot_std)
    if options.plot and bootstraps > 1:
        plot_fit(fn, cor, tmin, tmax, options, boot_averages, errors=boot_std)

    if return_chi:
        return boot_averages, boot_std, chi_sqr/dof
    if return_quality:
        return boot_averages, boot_std, quality_of_fit(dof, chi_sqr)
    else:
        return boot_averages, boot_std


def quality_of_fit(degrees_of_freedom, chi_sqr):
    dof = degrees_of_freedom
    return gammaincc(dof/2.0, chi_sqr / 2.0)

def write_fitted_cor(fn, cor, tmin, tmax, options, fitted_params, errors=None, postfix=None):
    if errors is None:
        fitted_errors = [0.0]* len(fitted_params)
    else:
        fitted_errors = errors
    if postfix:
        filestub = options.output_stub + postfix
    else:
        filestub = options.output_stub
    header="#fit {}, ({},{}), {}, {}, {}".format(fn.description, tmin, tmax, np.array(fitted_params), np.array(fitted_errors), options.period)
    cor.writeasv(filestub+".fittedcor.out", header=header)
    header="#fit_emass {}, ({},{}), {}, {}, {}".format(fn.description, tmin, tmax, np.array(fitted_params), np.array(fitted_errors), options.period)
    cor.writeemass(filestub+".fittedemass.out", dt=1, header=header)


def plot_fit(fn, cor, tmin, tmax, options, fitted_params, errors=None, postfix=None):
    import plot_helpers
    emass_dt = 1

    X = np.linspace(tmin, tmax, 200 * 5)

    if errors is None:
        fitted_errors = [0.0]* len(fitted_params)
    else:
        fitted_errors = errors

    fig = plt.figure()
    corplot = plt.subplot(211)
    cordata = corplot.errorbar(cor.times, cor.average_sub_vev().values(),
                               yerr=cor.jackknifed_errors().values(), fmt='o')
    corfit, = corplot.plot(X, fn.formula(fitted_params, X), lw=2.0)
    single_fit = None
    if fn.description == "exp" or "subtract" in fn.description:
        corplot.legend([cordata, corfit], ["Correlator data", fn.template.format(*fitted_params)], loc='best')
    else:
        single = functions["single_exp"]()
        single_fit, = corplot.plot(X, single.formula(fitted_params[:2], X), ls="-.", lw=2.0)
        corplot.legend([cordata, corfit, single_fit], ["Correlator data", fn.template.format(*fitted_params), "single_exp with these values"], loc='best')

    corplot.set_ylabel("Fit Correlator")

    corvals = cor.average_sub_vev().values()

    plt.ylim(plot_helpers.auto_fit_range(min(corvals),max(corvals)))
    plt.xlim([0, tmax + 2])
    emass = cor.periodic_effective_mass(emass_dt, fast=False, period=options.period)
    emass_errors = cor.periodic_effective_mass_errors(emass_dt, fast=False, period=options.period).values()
    emassplot = plt.subplot(212)
    emassplot.set_ylabel("${\mathrm{\mathbf{m}_{eff}}}$")
    dataplt = emassplot.errorbar(emass.keys(), emass.values(), yerr=emass_errors, fmt='o')
    named_params = {n: (m, e) for n, m, e in zip(fn.parameter_names, fitted_params, fitted_errors)}
    mass, mass_err = named_params["mass"]
    # abovefitline = emassplot.plot(range(tmin, tmax+1), [mass+mass_err]*len(range(tmin, tmax+1)), ls="dashed", color="b")
    fitplt, = emassplot.plot(range(tmin, tmax+1), [mass]*len(range(tmin, tmax+1)), ls="dotted", color="r")
    # belowfitline = emassplot.plot(range(tmin, tmax+1), [mass-mass_err]*len(range(tmin, tmax+1)), ls="dashed", color="b")
    fitpoints = fn.formula(fitted_params, np.arange(tmin, tmax+1))
    emassfit = []
    emassfit_range = []
    dt = emass_dt
    for i in range(len(fitpoints))[:-dt]:

        try:
            fitemass = (1.0 / float(dt)) * math.acosh((fitpoints[i+dt] + fitpoints[i-dt])/(2.0*fitpoints[i]))
            emassfit.append(fitemass)
            emassfit_range.append(tmin+i)
        except ValueError:
            fitemass = 0.0
    emass_fit = emassplot.plot(emassfit_range, emassfit, color="k")

    emassplot.legend([dataplt, fitplt], ["Emass of data", u"fit mass={:.5f}\xb1{:.5f}".format(mass, mass_err)], loc='best')
    # plt.ylim([min(min(emass.values()),-0.01), max(emass.values())*1.2])
    plt.xlim([0, tmax + 2])

    if options.output_stub:
        if postfix:
            filestub = options.output_stub + postfix
        else:
            filestub = options.output_stub
        logging.info("Saving plot to {}".format(filestub+".png"))
        fig.set_size_inches(18.5, 10.5)
        plt.savefig(filestub+".png")
    else:
        plt.show()
    plt.close()


def plot_histograms(names, paramters, options):
    logging.info("Plotting histograms of the bootstrap fits")
    import histo
    for index, name in enumerate(names):
        data = [p[index] for p in paramters]
        if options.output_stub:
            stub = "{}.{}.histo".format(options.output_stub, name)
        else:
            stub = None
        histo.make_histogram(data, options, stub, 100)



def bootstrap_cfgs(cor):
    return np.random.choice(cor.configs, size=len(cor.configs))


def bootstrap(cor, filelog=None):
    newcfgs = bootstrap_cfgs(cor)
    if filelog:
        with open(filelog+".straps", 'a') as bootfile:
            bootfile.write(",".join([str(c) for c in newcfgs]))
            bootfile.write("\n")
    newcor = {i: cor.get(config=c) for i, c in enumerate(newcfgs)}
    if cor.vev1 is None:
        newvev1 = newvev2 = None
    else:
        newvev1 = {i: cor.vev1[c] for i, c in enumerate(newcfgs)}
        newvev2 = {i: cor.vev2[c] for i, c in enumerate(newcfgs)}
    return correlator.Correlator.fromDataDicts(newcor, newvev1, newvev2)


def jackknife_ensemble(cor):
    oldcfgs = cor.configs
    jkensambles = []
    for cfg in cor.configs:
        newcfgs = [n for n in cor.configs if n != cfg]
        newcor = {i: cor.get(config=c) for i, c in enumerate(newcfgs)}
        if cor.vev1 is None:
            newvev1 = newvev2 = None
        else:
            newvev1 = {i: cor.vev1[c] for i, c in enumerate(newcfgs)}
            newvev2 = {i: cor.vev2[c] for i, c in enumerate(newcfgs)}
        # jkensambles.append()
        yield correlator.Correlator.fromDataDicts(newcor, newvev1, newvev2)

def bootstrap_ensamble(cor, N=NBOOTSTRAPS, filelog=None, jackknife=False):
    if jackknife:
        logging.warn("using jackknife instead of bootstrap")
        return jackknife_ensemble(cor)

    if N > 1:
        if filelog:
            with open(filelog+".straps", 'w') as bootfile:
                bootfile.write("# boot straps used for fitting")
        return (bootstrap(cor, filelog) for i in range(N))
    else:
        logging.info("Not bootstraping!")
        return [cor]


def covariance_matrix(cor, tmin, tmax):
    nm1 = (1.0 / (len(cor.configs) - 1))
    nm0 = 1.0 / (len(cor.configs))
    mymat = np.zeros(((tmax - tmin), (tmax - tmin)))
    start_time = cor.times[0]
    aoc = np.fromiter(cor.average_sub_vev().values(), np.float)[tmin-start_time:tmax-start_time]
    for v in cor.data.values():
        b = np.fromiter(v.values(), np.float)[tmin-start_time:tmax-start_time] - aoc
        # b = np.array(v.values()[tmin-start_time:tmax-start_time]).flat-aoc
        mymat += np.outer(b, b)
    return mymat*nm1*nm0


def jk_covariance_matrix(cor, tmin, tmax):
    """
    Here is the formula for the covariance matrix I am using:
    Let  < O >_J =  1/(N-1)  sum over bins, excluding Jth bin    O[bin]    for  N bins
    then   < C(t) >_J = < O(t) Obar(0) >_J - < O >_J < Obar >_J      (complex values, of course)
    Next, let  < C(t) >_avg = 1/N  *  sum over J < C(t) >_J
    then  Cov(t,t') = (1 - 1/N)  sum over J   ( <C(t)>_J - <C(t)>_avg ) *   ( <C(t')>_J - <C(t')>_avg )
    The last expression above is the jackknife estimate of the covariance.

    """
    mymat = np.zeros(((tmax - tmin), (tmax - tmin)))

    N = len(cor.configs)
    logging.debug("making covariance matrix with jackknifes")
    for t1 in range(tmin, tmax):
        for t2 in range(tmin, tmax):
            index1 = t1-tmin
            index2 = t2-tmin
            O_i = [v[t1] for k,v in cor.jackknife_average_sub_vev().iteritems()]
            O_j = [v[t2] for k,v in cor.jackknife_average_sub_vev().iteritems()]
            mean1 = np.mean(O_i)
            mean2 = np.mean(O_j)
            mymat[(index1,index2)] = (N-1) * np.mean( [(i-mean1)*(j-mean2)  for i,j in zip(O_i,O_j)])
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
        B[j][j] = 1.0 / tjj**2 - S / tjj  # noqa
        for i in reversed(range(j)):
            B[j][i] = B[i][j] = -np.sum([t[i][k] * B[k][j] for k in range(i + 1, nrows)]) / t[i][i]
    return B


def best_fit_range(fn, cor, options=None):
    logging.info("Finding best fit range")
    logging.debug("Temporarily setting the logger to warnings only")
    logger = logging.getLogger()
    previous_loglevel = logger.level
    logger.setLevel(ALWAYSINFO)
    best = 0
    best_ranges = []
    for tmin in cor.times:
        if fn.subtract and tmin == min(cor.times) or tmin < 1:
            continue
        tmaxes = [options.time_end] if options.time_end else range(tmin + len(fn.parameter_names)*4, max(cor.times))
        for tmax in tmaxes:
            if tmin > tmax:
                continue
            try:
                _, _, qual = fit(fn, cor, tmin, tmax, filestub=None, bootstraps=1, return_chi=True, return_quality=False, options=options)
                metric = 1/qual
                if qual < 1.5:
                    metric = (tmax-tmin)+metric
                # else:
                #     metric = qual
                #if metric > best:
                # best = metric
                best_ranges.append((metric, tmin, tmax))
                if metric > 0.2:
                    logging.log(ALWAYSINFO, "Fit range ({},{})"
                                " is good with chi/dof {} using {} points".format(tmin, tmax, qual, tmax-tmin))
            except RuntimeError:
                logging.warn("Fitter failed, skipping this tmin,tmax {},{}".format(tmin, tmax))
            except InversionError:
                logging.warn("Covariance matrix failed, skipping this tmin,tmax {},{}".format(tmin, tmax))
            # except Exception:
            #     logging.warn("Fitter failed, skipping this tmin,tmax")
    logger.setLevel(previous_loglevel)
    logging.debug("Restored logging state to original")
    return [(tmin, tmax) for _, tmin, tmax in sorted(best_ranges, reverse=True)]


def auto_fit(funct, cor, filestub=None, bootstraps=NBOOTSTRAPS, return_quality=False, options=None):
    fit_ranges = best_fit_range(funct, cor, options=options)
    logging.info("trying ranges {}".format(repr(fit_ranges)))
    for tmin, tmax in fit_ranges:
        logging.info("Trying fit range {}, {}".format(tmin, tmax))
        try:
            results = fit(funct, cor, tmin, tmax, filestub=filestub,
                          bootstraps=bootstraps, options=options)
            logging.info("Auto Fit sucessfully!")
            return  # (tmin, tmax) + results  # Need to return what fit range was done
        except RuntimeError as e:
            logging.warn("Fit range {} {} failed, trying next best".format(tmin, tmax))
            logging.warn("errored with {}".format(e))
            continue
        except InversionError as e:
            logging.warn("Fit range {} {} failed, trying next best".format(tmin, tmax))
            logging.warn("errored with {}".format(e))
        continue


def allfits(funct, cor, filestub=None, bootstraps=NBOOTSTRAPS, options=None):
    logging.info("Fitting ALL fit ranges")
    for tmin in cor.times:
        tmaxes = [options.time_end] if options.time_end else range(tmin + len(funct.parameter_names)*4, max(cor.times))
        for tmax in tmaxes:
            if tmin > tmax:
                continue
            try:
                _, _, qual = fit(funct, cor, tmin, tmax, filestub=filestub+"_{}_{}".format(tmin, tmax), bootstraps=bootstraps, return_chi=True, return_quality=False, options=options)
            except RuntimeError:
                logging.warn("Fitter failed, skipping this tmin,tmax {},{}".format(tmin, tmax))
            except InversionError:
                logging.warn("Covariance matrix failed, skipping this tmin,tmax {},{}".format(tmin, tmax))
    logging.info("fit all ranges")
    exit(0)


class InversionError(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


def bestInverse(M):
    TOLERANCE = 1.E-5

    def invert_error(i):
        return np.max(np.abs((np.dot(M, i) - np.identity(len(i)))))

    try:
        inv = linalg.inv(M)         # Will raise error if singular
    except np.linalg.linalg.LinAlgError:
        raise InversionError("Could not invert within tolerance")

    error = invert_error(inv)

    try:
        chol = linalg.cholesky(M, check_finite=False)
    except np.linalg.linalg.LinAlgError:
        logging.error("Not positive definite!")
        logging.exception("Could not invert Not positive definite!")
        raise InversionError("Could not invert within tolerance")

    else:
        chol_inv = CholeskyInverse(chol)
        chol_error = invert_error(chol_inv)

        logging.debug("Inversions errors were inv={}, chol={}".format(error, chol_error))

        if chol_error > TOLERANCE and error > TOLERANCE:
            logging.error("Error, {}".format(error))
            logging.error("chol_Error, {}".format(chol_error))
            raise InversionError("Could not invert within tolerance")

        if chol_error < error:
            logging.debug("Using choleskey inverse")
            inv = chol_inv
        else:
            logging.debug("Using standard inverse")
    return inv


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="compute fits", parents=[fitparser])

    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.DEBUG)
        logging.debug("Verbose debuging mode activated")
    else:
        logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.INFO)

    if args.output_stub is not None:
        root = logging.getLogger()
        errfilename = args.output_stub+".err"
        errfilehandler = logging.FileHandler(errfilename, delay=True)
        errfilehandler.setLevel(logging.WARNING)
        formatter = logging.Formatter('%(levelname)s: %(message)s')
        errfilehandler.setFormatter(formatter)
        root.addHandler(errfilehandler)
        logfilename = args.output_stub+".log"
        logfilehandler = logging.FileHandler(logfilename, delay=True)
        logfilehandler.setLevel(logging.INFO)
        formatter = logging.Formatter('%(levelname)s: %(message)s')
        logfilehandler.setFormatter(formatter)
        root.addHandler(logfilehandler)

    if args.output_stub:
        args.output_stub = os.path.splitext(args.output_stub)[0]
        outdir = os.path.dirname(args.output_stub)
        if not os.path.exists(outdir):
            logging.info("directory for output {} does not exist, atempting to create".format(outdir))
            if outdir is not "":
                os.makedirs(outdir)


    if args.output_stub and args.skip_done:
        filename = args.output_stub+".boot"
        try:
            if os.stat(filename).st_size > 0:
                logging.info(".boot file exists and not empty, skip fit")
                exit(0)
            else:
                logging.warn(".boot file exists but is empty!")
        except OSError as e:
            logging.info("running fit")


    if args.random:
        logging.info("Setting random seed to %s", args.random)
        np.random.seed(args.random)
#    print np.random.get_state()

    corrfile = args.inputfile
    vev1 = args.vev
    vev2 = vev1
    if args.vev2:
        vev2 = args.vev2


    cor = build_corr.corr_and_vev_from_pickle(corrfile, vev1, vev2)

    if args.symmetric or args.antisymmetric:
        corsym = cor.determine_symmetry()
        if corsym is None:
            logging.error("called with symmetric but correlator isnt")
            raise RuntimeError("called with symmetric but correlator isnt")
        logging.info("correlator found to be {}".format(corsym))
        cor.make_symmetric()




    if args.full:
        # check
        period = max(cor.times)+1
        if list(range(0, period)) != cor.times:
            raise RuntimeError("correlator times are not contiguous required by --full")
        args.period = period

    cor.prune_invalid(delete=True, sigma=args.prune)

    if args.bin:
        cor = cor.reduce_to_bins(args.bin)

    if not args.period:
        if cor.numconfigs == 551:
            logging.warning("period not set, guessing by confiigs, setting to 128")
            args.period = 128
        if cor.numconfigs == 412:
            logging.warning("period not set, guessing by confiigs, setting to 256")
            args.period = 256

    funct = functions[args.function](Nt=args.period)

    if args.alltimes:
        allfits(funct, cor, filestub=args.output_stub, bootstraps=args.bootstraps, options=args)

    if args.debugguess and (not args.time_start or not args.time_end):
        logging.warn("debugguess set without fit range, setting to max")
        args.maxrange = True

    if args.maxrange:
        logging.info("setting trange to the all of the data")
        args.time_start = min(cor.times)
        args.time_end = max(cor.times)
    if args.tmax:
        newtmax = max(cor.times)
        logging.info("setting tmax to {}".format(newtmax))
        args.time_end = newtmax

    tmin = args.time_start
    tmax = args.time_end
    fit_ranges = [(tmin, tmax)]
    if args.time_start is None:
        auto_fit(funct, cor, filestub=args.output_stub, bootstraps=args.bootstraps, options=args)
        exit()

    try:
        fit(funct, cor, tmin, tmax, filestub=args.output_stub, bootstraps=args.bootstraps, options=args)
    except InvalidFit:
        logging.error("Function does not have a fallback, fit failed")
        exit(-1)
        logging.error("Fit was invalid, trying backup")
        if funct.fallback and not args.nofallback:
            logging.error("function has a fallback {}".format(funct.fallback))
            fallback = functions[funct.fallback](Nt=args.period)
            fit(fallback, cor, tmin, tmax, filestub=args.output_stub, bootstraps=args.bootstraps, options=args)
        else:
            logging.error("Function does not have a fallback, fit failed")
            exit(-1)
