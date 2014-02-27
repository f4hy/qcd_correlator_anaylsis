#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
import logging
import correlator
import build_corr
import pylab
import argparse
import os

from fitfunctions import *  # noqa
import inspect
import sys

from scipy import linalg
from scipy import stats
from scipy.special import gammaincc
from scipy.optimize import leastsq
from scipy.optimize import fmin
from scipy.optimize import fmin_slsqp
from scipy.optimize import fmin_l_bfgs_b
# from scipy.optimize import minimize
OUTPUT = 25
ALWAYSINFO = 26
logging.addLevelName(OUTPUT, "OUTPUT")
logging.addLevelName(ALWAYSINFO, "INFO")

Nt = 128
NBOOTSTRAPS = 1000


def fit(fn, cor, tmin, tmax, filestub=None, bootstraps=NBOOTSTRAPS, return_quality=False, unsafe=False, return_chi=False, write_each_boot=None):
    if(tmax-tmin < len(fn.parameter_names)):
        raise ValueError("Can not fit to less points than parameters")

    results = logging.getLogger("results")
    if filestub and not results.handlers:
        filename = filestub+".log"
        filehandler = logging.FileHandler(filename)
        filehandler.level = OUTPUT
        results.addHandler(filehandler)
        logging.info("Writing output to file {}".format(filename))

    results.info("Fitting data to {} from t={} to t={} using {} bootstrap samples".format(
        fn.description, tmin, tmax, bootstraps))

    tmax = tmax+1  # I use ranges, so this needs to be offset by one
    fun = lambda v, mx, my: (fn.formula(v, mx) - my)

    initial_guess = fn.starting_guess(cor, tmax-1, tmin)
    logging.info("Starting with initial_guess: {}".format(repr(initial_guess)))
    x = np.array(range(tmin, tmax))
    orig_ave_cor = cor.average_sub_vev()
    y = [orig_ave_cor[t] for t in range(tmin, tmax)]
    original_ensamble_params, success = leastsq(fun, initial_guess, args=(x, y), maxfev=10000)
    if not success:
        raise ValueError()
    initial_guess = original_ensamble_params

    def cov_fit(correlator, guess):
        ave_cor = correlator.average_sub_vev()
        y = [ave_cor[t] for t in range(tmin, tmax)]
        cov = covariance_matrix(correlator, tmin, tmax)
        inv_cov = bestInverse(cov)
        start_time = correlator.times[0]
        aoc = np.fromiter(ave_cor.values(), np.float)[tmin-start_time:tmax-start_time]

        def cov_fun(g):
            """ Function to be minizied. computed using matrix mult"""
            vect = aoc - fn.formula(g, x)
            return vect.dot(inv_cov).dot(vect)
        uncorrelated_fit_values, success = leastsq(fun, guess, args=(x, y), maxfev=100000)
        if not success:
            raise ValueError()

        guess = uncorrelated_fit_values
        #results = fmin(cov_fun, fn.starting_guess, ftol=1.E-7, maxfun=1000000, maxiter=1000000, full_output=1, disp=0, retall=0)
        if guess[0] < 0.0:
            logging.warn("first pass fit value found mass to be negative {}, lets flip it".format(guess[0]))
            guess[0] = -guess[0]

        if len(guess) > 2 and guess[2] < 0.0:
            logging.warn("first pass fit value found mass2 to be negative {}, lets flip it".format(guess[2]))
            logging.info("first pass results are {}".format(repr(guess)))
            guess[2] = -guess[2]
        def clamp(n, minn, maxn):
                return max(min(maxn, n), minn)
        bounded_guess = [clamp(g,b[0],b[1]) for g,b in zip(guess,fn.bounds)]
        logging.debug("guess {}, bounded guess {}".format(repr(guess),repr(bounded_guess)))
        newresults = fmin_slsqp(cov_fun, bounded_guess, bounds=fn.bounds, full_output=1, disp=0, iter=10000)
        logging.debug("fit value {}".format(repr(newresults)))
        results = newresults
        covariant_fit, fit_info, flag = results[0], results[1:3], results[3]
        if covariant_fit[0] < 0.0:
            logging.error("Fitter gave negative mass {}!!! Error!".format(covariant_fit[0]))
            raise RuntimeError("Fitter sanity failed")

        #would like to use minimize, but seems to be not installed
        # covariant_fit = minimize(cov_fun, initial_guess)
        #logging.debug("Fit results: f() ={}, Iterations={}, Function evaluations={}".format(*fit_info))
        logging.debug("Fit results: f() ={}, Iterations={}".format(*fit_info))
        if flag != 0:
            logging.error("Fitter flag set to {}. Error!".format(flag))
            raise RuntimeError("Fitter failed")

        logging.debug("Covariant fit {}, regular fit {}".format(repr(covariant_fit),
                                                                repr(uncorrelated_fit_values)))

        return(covariant_fit)

    boot_params = []
    for strap in bootstrap_ensamble(cor, N=bootstraps, filelog=write_each_boot):
        newguess = fn.starting_guess(strap, tmax-1, tmin)
        fitted_params = cov_fit(strap, newguess)
        boot_params.append(fitted_params)

    original_ensamble_correlatedfit = cov_fit(cor, original_ensamble_params)

    results.info('')
    results.info('Uncorelated total fit: %s', {n: p for n, p in zip(fn.parameter_names, original_ensamble_params)})
    results.info('Correlated total fit:  %s', {n: p for n, p in zip(fn.parameter_names, original_ensamble_correlatedfit)})
    boot_averages = np.mean(boot_params, 0)
    boot_std = np.std(boot_params, 0)
    boota = np.array(boot_params)
    upper_quartiles = [stats.scoreatpercentile(boota[:,i],75) for i in range(len(boot_averages))]
    medians = [stats.scoreatpercentile(boota[:,i],50) for i in range(len(boot_averages))]
    lower_quartiles = [stats.scoreatpercentile(boota[:,i],25) for i in range(len(boot_averages))]
    inter_range = [stats.scoreatpercentile(boota[:,i],75) - stats.scoreatpercentile(boota[:,i],25) for i in range(len(boot_averages))]


    for name, boot, original, err in zip(fn.parameter_names, boot_averages, original_ensamble_correlatedfit, boot_std):
        bias = abs(boot-original)
        percent_bias = abs(boot-original)/original
        results.info('Bootstrap Bias in {:<10}: {:.3%}'.format(name, percent_bias))
        if bias > err*2:
            results.error('Bootstrap Bias in {:<10}: {:.3%}'.format(name, percent_bias))
            results.error("Bootstrap average does not agree with ensamble average!"
                          "\nNot enough statistics for this for to be valid!!!\n")
            if not unsafe:
                results.critical("Exiting! Run with --unsafe to fit anyway")
                raise RuntimeError("Bootstrap average does not agree with ensamble average")

    for name, ave, med, std, iqr in zip(fn.parameter_names, boot_averages, medians, boot_std, inter_range):
        skew = abs(ave-med)/ave
        dist_skew = abs(std-iqr)/iqr
        if skew > 1.0:
            results.error("for {} diff of bootstrap average and bootstrap median is {:.3%}".format(name, skew))
            results.error("Bootstrap distrubtion is skewed!!")
        else:
            results.info("for {} diff of bootstrap average and bootstrap median is {:.3%}".format(name, skew))
        if dist_skew > 1.0:
            results.error("for {} diff of standard deviation and interquartile range is {:.3%}".format(name, dist_skew))
            results.error("Large outliers present in bootstrap fits!!")
        else:
            results.info("for {} diff of standard deviation and interquartile range is {:.3%}".format(name, dist_skew))


    results.info("")
    results.log(OUTPUT, "Full esamble fitted parameters t=%2d to %2d---------------------", tmin, tmax)
    results.log(OUTPUT, "Name      : Average,        STD,           (1st Quart, Median, 3rd Quart, IQR)")
    for name, ave, std, low, med, up, iqr in zip(fn.parameter_names, original_ensamble_correlatedfit, boot_std,
                                            upper_quartiles, medians, lower_quartiles, inter_range):
        results.log(OUTPUT, u"{:<10}: {:<15.10f} \u00b1 {:<10g}   ({:<9.6f}, {:<9.6f}, {:<9.6f}, {:<9.6f})".format(name, ave, std, low, med, up, iqr))
    results.log(OUTPUT, "--------------------------------------------------------")

    v = original_ensamble_correlatedfit
    cov = covariance_matrix(cor, tmin, tmax)
    inv_cov = bestInverse(cov)
    chi_sqr = np.sum(((orig_ave_cor[t] - fn.formula(v, t)) * inv_cov[t - tmin][tp - tmin] * (orig_ave_cor[tp] - fn.formula(v, tp))
                      for t in range(tmin, tmax) for tp in range(tmin, tmax)))

    dof = len(x) - len(fn.parameter_names)
    results.log(OUTPUT, u'\u03c7\u00b2 ={},   \u03c7\u00b2 / dof = {}, Qual {}\n'.format(
        chi_sqr, chi_sqr/dof, quality_of_fit(dof, chi_sqr)))

    if write_each_boot:
        results.info("writing each bootstrap to {}.boot".format(write_each_boot))
        with open(write_each_boot+".boot", 'w') as bootfile:
            str_ensamble_params = ", ".join([str(p) for p in original_ensamble_params])
            bootfile.write("#bootstrap, {}, \t ensamble mean: {}\n".format(", ".join(fn.parameter_names), str_ensamble_params))
            for i, params in enumerate(boot_params):
                strparams = ", ".join([str(p) for p in params])
                bootfile.write("{}, {}\n".format(i, strparams))

    if return_chi:
        return original_ensamble_correlatedfit , boot_std, chi_sqr/dof
    if return_quality:
        return original_ensamble_correlatedfit, boot_std, quality_of_fit(dof, chi_sqr)
    else:
        return original_ensamble_correlatedfit, boot_std


def quality_of_fit(degrees_of_freedom, chi_sqr):
    dof = degrees_of_freedom
    return gammaincc(dof/2.0, chi_sqr / 2.0)


def plot_fit(fn, cor, tmin, tmax, filestub=None, bootstraps=NBOOTSTRAPS, unsafe=False):
    emass_dt = 3

    X = np.linspace(tmin, tmax, 200 * 5)
    massX = np.linspace(tmin, tmax-emass_dt, 200 * 5)
    fitted_params, fitted_errors = fit(fn, cor, tmin, tmax, filestub, bootstraps=bootstraps, unsafe=unsafe)

    plt.figure()
    corplot = plt.subplot(211)
    cordata = corplot.errorbar(cor.times, cor.average_sub_vev().values(),
                               yerr=cor.jackknifed_errors().values(), fmt='o')
    corfit, = corplot.plot(X, fn.formula(fitted_params, X))
    corplot.legend([cordata, corfit], ["data", fn.template.format(*fitted_params)])
    plt.ylim([0, max(cor.average_sub_vev().values())])
    plt.xlim([0, tmax + 2])
    emass = cor.effective_mass(emass_dt)
    emass_errors = cor.effective_mass_errors(emass_dt).values()
    emassplot = plt.subplot(212)
    dataplt = emassplot.errorbar(emass.keys(), emass.values(), yerr=emass_errors, fmt='o')
    named_params = {n: (m, e) for n, m, e in zip(fn.parameter_names, fitted_params, fitted_errors)}
    mass, mass_err = named_params["mass"]
    fitplt = emassplot.errorbar(massX, mass * np.ones_like(massX), yerr=mass_err)
    plt.legend([dataplt, fitplt], ["data", u"fit mass={:.5f}\xb1{:.5f}".format(mass, mass_err)])
    plt.ylim([0, max(emass.values())*1.2])
    plt.xlim([0, tmax + 2])
    if(filestub):
        logging.info("Saving plot to {}".format(filestub+".png"))
        plt.savefig(filestub+".png")
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


def best_fit_range(fn, cor):
    logging.info("Finding best fit range")
    logging.debug("Temporarily setting the logger to warnings only")
    logger = logging.getLogger()
    previous_loglevel = logger.level
    logger.setLevel(ALWAYSINFO)
    best = 100
    best_ranges = []
    for tmin in cor.times:
        for tmax in range(tmin + 4, max(cor.times)):
            try:
                _, _, chi = fit(fn, cor, tmin, tmax, filestub=None, bootstraps=1, return_chi=True)
                metric = abs(chi-1.0)
                # if metric < best:
                best = metric
                best_ranges.append((metric, tmin, tmax))
                if best < 1.0:
                    logging.log(ALWAYSINFO, "Fit range ({},{})"
                                " is good with chi/dof {}".format(tmin, tmax, chi))
            except RuntimeError:
                logging.warn("Fitter failed, skipping this tmin,tmax {},{}".format(tmin,tmax))
            # except Exception:
            #     logging.warn("Fitter failed, skipping this tmin,tmax")
    logger.setLevel(previous_loglevel)
    logging.debug("Restored logging state to original")
    return [(tmin, tmax) for _, tmin, tmax in sorted(best_ranges)]


def auto_fit(funct, cor, filestub=None, bootstraps=NBOOTSTRAPS, return_quality=False, unsafe=False, write_each_boot=None):
    fit_ranges = best_fit_range(funct, cor)
    for tmin, tmax in fit_ranges:
        logging.info("Trying fit range {}, {}".format(tmin, tmax))
        try:
            results = fit(funct, cor, tmin, tmax, filestub=filestub,
                          bootstraps=bootstraps, unsafe=unsafe, return_quality=return_quality, write_each_boot=write_each_boot)
            logging.info("Auto Fit sucessfully!")
            return (tmin, tmax) + results  # Need to return what fit range was done
        except RuntimeError:
            logging.warn("Fit range {} {} failed, trying next best".format(tmin, tmax))
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

args = None


if __name__ == "__main__":

    function_list = inspect.getmembers(sys.modules["fitfunctions"], inspect.isclass)
    functions = {name: f for name, f in function_list}

    parser = argparse.ArgumentParser(description="compute fits")
    parser.add_argument("-i", "--inputfile", type=str, required=True,
                        help="Correlator file to read from")
    parser.add_argument("-o", "--output_stub", type=str, required=False,
                        help="stub of name to write output to")
    parser.add_argument("-wb", "--write_boot", type=str, required=False,
                        help="stub of name to write each bootstrap output to")
    parser.add_argument("-v1", "--vev", type=str, required=False,
                        help="vev file to read from")
    parser.add_argument("-v2", "--vev2", type=str, required=False,
                        help="vev2 file to read from")
    parser.add_argument("-ts", "--time-start", type=int, required=False,
                        help="first time slice to start a fit, can be a list of times")
    parser.add_argument("-te", "--time-end", type=int, required=False,
                        help="last time slice to fit, can be a list of times")
    parser.add_argument("-b", "--bootstraps", type=int, required=False, default=NBOOTSTRAPS,
                        help="Number of straps")
    parser.add_argument("-p", "--plot", action="store_true", required=False,
                        help="Plot the resulting fit")
    parser.add_argument("-Nt", "--period", type=int, required=False,
                        help="Period in time direction (not required for all functions)")
    parser.add_argument("-r", "--random", type=int, required=False,
                        help="set the random seed")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="increase output verbosity")
    parser.add_argument("--unsafe", action="store_true",
                        help="Code usually exits if something goes wrong "
                        "This option will cause the code to fit anyway.")
    parser.add_argument("-f", "--function", choices=functions.keys(),
                        required=False, default="periodic_exp", help="function to fit to")

    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.DEBUG)
        logging.debug("Verbose debuging mode activated")
    else:
        logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.INFO)

    funct = functions[args.function](Nt=args.period)

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
    cor.prune_invalid(delete=True)
    tmin = args.time_start
    tmax = args.time_end
    fit_ranges = [(tmin, tmax)]
    if not args.time_start:
        print args.output_stub
        auto_fit(funct, cor, filestub=args.output_stub, bootstraps=args.bootstraps, unsafe=args.unsafe, write_each_boot=args.write_boot)
        exit()

    if args.plot:
        plot_fit(funct, cor, tmin, tmax, filestub=args.output_stub,
                 bootstraps=args.bootstraps, unsafe=args.unsafe)
    else:
        fit(funct, cor, tmin, tmax, filestub=args.output_stub,
            bootstraps=args.bootstraps, unsafe=args.unsafe,
            write_each_boot=args.write_boot)
