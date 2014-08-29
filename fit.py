#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
import logging
import correlator
import build_corr
import argparse
import os

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
OUTPUT = 25
ALWAYSINFO = 26
logging.addLevelName(OUTPUT, "OUTPUT")
logging.addLevelName(ALWAYSINFO, "INFO")

Nt = 128
NBOOTSTRAPS = 1000


def fit(fn, cor, tmin, tmax, filestub=None, bootstraps=NBOOTSTRAPS, return_quality=False,
        return_chi=False, options=None):
    if(tmax-tmin < len(fn.parameter_names)):
        raise InvalidFit("Can not fit to less points than parameters")

    results = logging.getLogger("results")
    if filestub and not results.handlers:
        filename = filestub+".log"
        filehandler = logging.FileHandler(filename)
        filehandler.level = OUTPUT
        results.addHandler(filehandler)
        logging.info("Writing output to file {}".format(filename))

    results.info("Fitting data to {} from t={} to t={} using {} bootstrap samples".format(
        fn.description, tmin, tmax, bootstraps))

    #tmax = tmax+1  # I use ranges, so this needs to be offset by one
    fitrange = range(tmin, tmax+1)
    fun = lambda v, mx, my: (fn.formula(v, mx) - my)
    initial_guess = fn.starting_guess(cor, tmax, tmin)
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
        return initial_guess, [0.01, 0.01, 0.01, 0.01]  # For testing initila guess in plot
    if not success:
        raise InvalidFit("original exnamble leastsq failed")
    if options.first_pass:
        initial_guess = original_ensamble_params
        logging.info("initial_guess after first pass: {}".format(repr(initial_guess)))

    def cov_fit(correlator, guess):
        ave_cor = correlator.average_sub_vev()
        y = [ave_cor[t] for t in fitrange]
        cov = covariance_matrix(correlator, tmin, tmax+1)
        inv_cov = bestInverse(cov)
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
        if m.get_fmin().is_valid:
            return minuit_results
        else:
            logging.error("minuit failed!!")
            logging.error("was at {}".format(minuit_results))
            return None

    # end cov_fit

    original_ensamble_correlatedfit = cov_fit(cor, initial_guess)
    isvalidfit = fn.valid(original_ensamble_correlatedfit)
    if not isvalidfit:
        raise InvalidFit("Full ensamble failed")


    boot_params = []
    failcount = 0
    attempted = 0
    for strap in bootstrap_ensamble(cor, N=bootstraps, filelog=filestub):
        attempted +=1
        if options.reguess:
            newguess = fn.starting_guess(strap, tmax, tmin)
        else:
            newguess = initial_guess
        fitted_params = cov_fit(strap, newguess)
        if fitted_params is not None:
            boot_params.append(fitted_params)
            logging.debug("bootstrap converged")
        else:
            logging.error("bootstrap failed to converge!")
            raise InvalidFit("one bootstrap failed")
            #raw_input("test")
            failcount+=1
            logging.debug("fails:{} attempts:{}, ratio:{}".format(failcount, attempted, failcount/float(attempted)))
            if failcount/float(attempted) > 0.15 and attempted > 40:
                raise InvalidFit("more than 20% of boostraps failed to converge")

    if failcount > 0:
        logging.warn("{} bootstraps did not converge!".format(bootstraps-len(boot_params)))
    if len(boot_params) < bootstraps * 0.9:
        logging.error("More that 10% of the straps failed")
        raise InvalidFit("more than 10% of boostraps failed to converge")

    results.info('')
    results.info('Uncorelated total fit: %s', {n: p for n, p in zip(fn.parameter_names, original_ensamble_params)})
    results.info('Correlated total fit:  %s', {n: p for n, p in zip(fn.parameter_names, original_ensamble_correlatedfit)})
    boot_averages = np.mean(boot_params, 0)
    boot_std = np.std(boot_params, 0)
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
        if dist_skew > 1.0:
            results.error("for {} diff of stddev and IQR is {:.3%}".format(name, dist_skew))
            results.error("Large outliers present in bootstrap fits!!")
            if not options.unsafe:
                results.critical("Exiting! Run with --unsafe to fit anyway")
                raise InvalidFit("Bootstrap average does not agree with ensamble average")
        else:
            results.info("for {} diff of standard deviation"
                         " and interquartile range is {:.3%}".format(name, dist_skew))

    results.info("")
    results.log(OUTPUT, "Full esamble fitted parameters t=%2d to %2d---------------------", tmin, tmax)
    results.log(OUTPUT, "Name      : Average,        STD,           (1st Quart, Median, 3rd Quart, IQR)")
    for name, ave, std, low, med, up, iqr in zip(fn.parameter_names, original_ensamble_correlatedfit, boot_std,
                                                 upper_quartiles, medians, lower_quartiles, inter_range):
        results.log(OUTPUT, u"{:<10}: {:<15.10f} \u00b1 {:<10g}   ({:<9.6f}, {:<9.6f}, {:<9.6f}, {:<9.6f})".format(name, ave, std, low, med, up, iqr))
    results.log(OUTPUT, "--------------------------------------------------------")

    v = original_ensamble_correlatedfit
    cov = covariance_matrix(cor, tmin, tmax+1)
    inv_cov = bestInverse(cov)
    chi_sqr = np.sum(((orig_ave_cor[t] - fn.formula(v, t)) * inv_cov[t - tmin][tp - tmin] * (orig_ave_cor[tp] - fn.formula(v, tp))
                      for t in fitrange for tp in fitrange))

    dof = len(x) - len(fn.parameter_names)
    results.log(OUTPUT, u'\u03c7\u00b2 ={},   \u03c7\u00b2 / dof = {}, Qual {}\n'.format(
        chi_sqr, chi_sqr/dof, quality_of_fit(dof, chi_sqr)))

    valid = True
    if original_ensamble_correlatedfit[0]*2 < min(cor.effective_mass(3).values()):
        logging.error("fit value much lower than effective mass {} , {}!!".format(original_ensamble_correlatedfit[0],
                                                                                  min(cor.effective_mass(3).values())))
        valid = False

    if options.write_each_boot and valid and filestub:
        results.info("writing each bootstrap to {}.boot".format(filestub))
        with open(filestub+".boot", 'w') as bootfile:
            str_ensamble_params = ", ".join([str(p) for p in original_ensamble_params])
            bootfile.write("#bootstrap, {}, \t ensamble mean: {}\n".format(", ".join(fn.parameter_names), str_ensamble_params))
            for i, params in enumerate(boot_params):
                strparams = ", ".join([str(p) for p in params])
                bootfile.write("{}, {}\n".format(i, strparams))


    if return_chi:
        return original_ensamble_correlatedfit, boot_std, chi_sqr/dof
    if return_quality:
        return original_ensamble_correlatedfit, boot_std, quality_of_fit(dof, chi_sqr)
    else:
        return original_ensamble_correlatedfit, boot_std


def quality_of_fit(degrees_of_freedom, chi_sqr):
    dof = degrees_of_freedom
    return gammaincc(dof/2.0, chi_sqr / 2.0)


def plot_fit(fn, cor, tmin, tmax, filestub=None, bootstraps=NBOOTSTRAPS, options=None):
    emass_dt = 3

    X = np.linspace(tmin, tmax, 200 * 5)
    fitted_params, fitted_errors = fit(fn, cor, tmin, tmax, filestub, bootstraps=bootstraps, options=options)

    fig = plt.figure()
    corplot = plt.subplot(211)
    cordata = corplot.errorbar(cor.times, cor.average_sub_vev().values(),
                               yerr=cor.jackknifed_errors().values(), fmt='o')
    corfit, = corplot.plot(X, fn.formula(fitted_params, X), lw=2.0)
    single_fit = None
    if fn.description == "exp" or "subtract" in fn.description:
        corplot.legend([cordata, corfit], ["Correlator data", fn.template.format(*fitted_params)])
    else:
        single = functions["single_exp"]()
        single_fit, = corplot.plot(X, single.formula(fitted_params[:2], X), ls="-.", lw=2.0)
        corplot.legend([cordata, corfit, single_fit], ["Correlator data", fn.template.format(*fitted_params), "single_exp with these values"])

    corplot.set_ylabel("Fit Correlator")

    corvals = cor.average_sub_vev().values()
    plt.ylim([min(min(corvals), 0), max(corvals)])
    plt.xlim([0, tmax + 2])
    emass = cor.effective_mass(emass_dt)
    emass_errors = cor.effective_mass_errors(emass_dt).values()
    emassplot = plt.subplot(212)
    emassplot.set_ylabel("${\mathrm{\mathbf{m}_{eff}}}$")
    dataplt = emassplot.errorbar(emass.keys(), emass.values(), yerr=emass_errors, fmt='o')
    named_params = {n: (m, e) for n, m, e in zip(fn.parameter_names, fitted_params, fitted_errors)}
    mass, mass_err = named_params["mass"]
    # abovefitline = emassplot.plot(range(tmin, tmax+1), [mass+mass_err]*len(range(tmin, tmax+1)), ls="dashed", color="b")
    fitplt = emassplot.plot(range(tmin, tmax+1), [mass]*len(range(tmin, tmax+1)), ls="dotted", color="r")
    # belowfitline = emassplot.plot(range(tmin, tmax+1), [mass-mass_err]*len(range(tmin, tmax+1)), ls="dashed", color="b")
    fitpoints = fn.formula(fitted_params, np.arange(tmin, tmax+1))
    emassfit = []
    dt = emass_dt
    for i in range(len(fitpoints))[:-dt]:
        fitemass = (1.0 / float(dt)) * np.log(fitpoints[i] / fitpoints[i + dt])
        emassfit.append(fitemass)
    emass_fit = emassplot.plot(range(tmin, tmax+1)[:-dt], emassfit)

    plt.legend([dataplt, fitplt], ["Emass of data", u"fit mass={:.5f}\xb1{:.5f}".format(mass, mass_err)])
    plt.ylim([min(min(emass.values()),0), max(emass.values())*1.2])
    plt.xlim([0, tmax + 2])

    if(filestub):
        logging.info("Saving plot to {}".format(filestub+".png"))
        fig.set_size_inches(18.5, 10.5)
        plt.savefig(filestub+".png")
        header="#fit {}, ({},{}), {}, {}".format(fn.description, tmin, tmax, np.array(fitted_params), np.array(fitted_errors))
        cor.writeasv(filestub+".fittedcor.out", header=header)
        header="#fit_emass {}, ({},{}), {}, {}".format(fn.description, tmin, tmax, np.array(fitted_params), np.array(fitted_errors))
        cor.writeemass(filestub+".fittedemass.out", header=header)
    else:
        plt.show()


def bootstrap_cfgs(cor):
    return np.random.choice(cor.configs, size=len(cor.configs))


def bootstrap(cor, filelog=None):
    newcfgs = bootstrap_cfgs(cor)
    if filelog:
        with open(filelog+".straps", 'a') as bootfile:
            bootfile.write(",".join([str(c) for c in newcfgs]))
            bootfile.write("\n")
    newcor = {i: cor.get(config=c) for i, c in enumerate(newcfgs)}
    newvev1 = {i: cor.vev1[c] for i, c in enumerate(newcfgs)}
    newvev2 = {i: cor.vev2[c] for i, c in enumerate(newcfgs)}
    return correlator.Correlator.fromDataDicts(newcor, newvev1, newvev2)


def bootstrap_ensamble(cor, N=NBOOTSTRAPS, filelog=None):
    if N > 1:
        with open(filelog+".straps", 'w') as bootfile:
            bootfile.write("# boot straps used for fitting")
        return [bootstrap(cor, filelog) for i in range(N)]
    else:
        logging.info("Not bootstraping!")
        return [cor]


def covariance_matrix(cor, tmin, tmax):
    nm1 = (1.0 / (len(cor.configs) - 1))
    nm0 = 1.0 / (len(cor.configs))
    mymat = np.zeros(((tmax - tmin), (tmax - tmin)))
    start_time = cor.times[0]
    aoc = np.fromiter(cor.average_over_configs().values(), np.float)[tmin-start_time:tmax-start_time]
    for v in cor.data.values():
        b = np.fromiter(v.values(), np.float)[tmin-start_time:tmax-start_time] - aoc
        # b = np.array(v.values()[tmin-start_time:tmax-start_time]).flat-aoc
        mymat += np.outer(b, b)
    return mymat*nm1*nm0


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
        if fn.subtract and tmin == min(cor.times):
            continue
        tmaxes = [options.time_end] if options.time_end else range(tmin + 4, max(cor.times))
        for tmax in tmaxes:
            try:
                _, _, qual = fit(fn, cor, tmin, tmax, filestub=None, bootstraps=1, return_chi=False, return_quality=True, options=options)
                #metric = abs(chi-1.0)
                metric = qual
                #if metric > best:
                # best = metric
                best_ranges.append((metric, tmin, tmax))
                if metric > 0.2:
                    logging.log(ALWAYSINFO, "Fit range ({},{})"
                                " is good with chi/dof {}".format(tmin, tmax, qual))
            except RuntimeError:
                logging.warn("Fitter failed, skipping this tmin,tmax {},{}".format(tmin, tmax))
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
            if args.plot:
                results = plot_fit(funct, cor, tmin, tmax, filestub=filestub,
                                   bootstraps=bootstraps, options=options)
            else:
                results = fit(funct, cor, tmin, tmax, filestub=filestub,
                              bootstraps=bootstraps, options=options)
            logging.info("Auto Fit sucessfully!")
            return  # (tmin, tmax) + results  # Need to return what fit range was done
        except RuntimeError as e:
            logging.warn("Fit range {} {} failed, trying next best".format(tmin, tmax))
            logging.warn("errored with {}".format(e))
            continue


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
        chol = linalg.cholesky(M, check_finite=False)
    except np.linalg.linalg.LinAlgError:
        logging.error("Not positive definite!")
        logging.exception("Could not invert Not positive definite!")
        raise

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

    parser = argparse.ArgumentParser(description="compute fits", parents=[fitparser])

    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.DEBUG)
        logging.debug("Verbose debuging mode activated")
    else:
        logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.INFO)

    if args.random:
        logging.info("Setting random seed to %s", args.random)
        np.random.seed(args.random)
#    print np.random.get_state()

    corrfile = args.inputfile
    vev1 = args.vev
    vev2 = vev1
    if args.vev2:
        vev2 = args.vev2

    if args.output_stub:
        args.output_stub = os.path.splitext(args.output_stub)[0]

    cor = build_corr.corr_and_vev_from_files_pandas(corrfile, vev1, vev2)
    cor.prune_invalid(delete=True, sigma=args.prune)

    if not args.period:
        if cor.numconfigs == 551:
            logging.warning("period not set, guessing by confiigs, setting to 128")
            args.period = 128
        if cor.numconfigs == 412:
            logging.warning("period not set, guessing by confiigs, setting to 256")
            args.period = 256

    funct = functions[args.function](Nt=args.period)

    if args.maxrange:
        logging.info("setting trange to the all of the data")
        args.time_start = min(cor.times)
        args.time_end = max(cor.times)
    if args.tmax:
        logging.info("setting tmax to {}".format(max(cor.times)))
        args.time_end = max(cor.times)
        
    tmin = args.time_start
    tmax = args.time_end
    fit_ranges = [(tmin, tmax)]
    if not args.time_start:
        print args.output_stub
        auto_fit(funct, cor, filestub=args.output_stub, bootstraps=args.bootstraps, options=args)
        exit()

    try:
        if args.plot:
            plot_fit(funct, cor, tmin, tmax, filestub=args.output_stub,
                     bootstraps=args.bootstraps, options=args)
        else:
            fit(funct, cor, tmin, tmax, filestub=args.output_stub, bootstraps=args.bootstraps, options=args)
    except InvalidFit:
        logging.error("Fit was invalid, trying backup")
        if funct.fallback and not args.nofallback:
            logging.error("function has a fallback {}".format(funct.fallback))
            fallback = functions[funct.fallback](Nt=args.period)
            if args.plot:
                plot_fit(fallback, cor, tmin, tmax, filestub=args.output_stub,
                         bootstraps=args.bootstraps, options=args)
            else:
                fit(fallback, cor, tmin, tmax, filestub=args.output_stub, bootstraps=args.bootstraps, options=args)
        else:
            logging.error("Function does not have a fallback, fit failed")
