#!/usr/bin/env python
import fit
import inspect
import sys
from simul_fitfunctions import *  # noqa
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




def consistant_argument_counts(args):
    number_of_files = len(args.inputfile)
    logging.info("Checking to make sure arguments are consistant with {} count".format(number_of_files))
    if args.vev:
        assert(number_of_files == len(args.vev))
        assert(number_of_files == len(args.vev2))
    else:
        args.vev=[None]*number_of_files
    if args.time_start:
        assert(number_of_files == len(args.time_start))
    if args.time_end:
        assert(number_of_files == len(args.time_end))


def mergecors(cors, tmins, tmaxs):
    configs = cors[0].configs
    assert(all(configs == c.configs for c in cors))


    newdata = {}
    dummyvev = {}
    boundries = set()
    for conf in configs:
        dummyvev[conf] = 0.0
        dataindex = 0
        confdata = {}
        for i in range(0,len(cors)):
            for t in range(tmins[i], tmaxs[i]+1):
                confdata[dataindex] = cors[i].get(config=conf, time=t)
                dataindex +=1
            boundries.add(dataindex-1)
            boundries.add(dataindex)
        newdata[conf] = confdata

    if cors[0].vev1 is None and cors[1].vev1 is None:
        dummyvev = None

    c = correlator.Correlator.fromDataDicts(newdata, dummyvev, dummyvev)
    c.emass_skip_times = boundries
    return c


def auto_fit(cors, options=None):
    logging.info("Finding best fit range")
    logging.debug("Temporarily setting the logger to warnings only")

    individual_fitfun = functions[options.function].individual(options.period)

    print cors
    ranges = []
    for cor in cors:
        ranges.append(fit.best_fit_range(individual_fitfun, cor, options)[0])

    print "best ranges"
    print ranges
    #exit()
    # logger = logging.getLogger()
    # previous_loglevel = logger.level
    # ALWAYSINFO = 26
    # logger.setLevel(ALWAYSINFO)
    # import itertools
    # tend = options.period/2
    # print tend

    # tends = [tend] * len(cors)



    # allowed_tmins = range(5,tend-5)
    # tmins = itertools.product(allowed_tmins, repeat=2)
    # tmins = [(x,y) for x,y in tmins if x>y]
    # best_ranges = []
    # for tmin in tmins:
    #     multicor = mergecors(cors, tmin, tends)
    #     funct = functions[options.function](Nt=options.period, ranges=zip(tmin, tends))
    #     try:
    #         _, _, qual = fit.fit(funct, multicor, min(multicor.times), max(multicor.times),
    #                              bootstraps=1, return_chi=False, return_quality=True, options=options)
    #         metric = qual
    #         best_ranges.append((metric,tmin))
    #     except RuntimeError:
    #         logging.warn("Fitter failed, skipping this tmin,tmax {},{}".format(*tmin))
    #     except fit.InversionError:
    #         logging.warn("Covariance matrix failed, skipping this tmin,tmax {},{}".format(*tmin))
    #     except Exception as e:
    #         logging.warn("Fitter failed error {}".format(e))
    #         logging.warn("Fitter failed, skipping this tmin,tmax {},{}".format(*tmin))

    # logger.setLevel(previous_loglevel)
    # logging.debug("Restored logging state to original")

    # sortedranges = [(metric, tmins) for metric, tmins in sorted(best_ranges, reverse=True)]
    # print sortedranges[:5]
    # exit()
    # for _, tmin in sortedranges[:10]:
    #     try:
    #         multicor = mergecors(cors, tmin, tends)
    #         funct = functions[options.function](Nt=options.period, ranges=zip(tmin, tends))
    #         fit.fit(funct, multicor, min(multicor.times), max(multicor.times),
    #                 filestub=options.output_stub, return_chi=False, return_quality=True, options=options)
    #     except RuntimeError:
    #         logging.warn("Fitter failed, skipping this tmin,tmax {},{}".format(*tmin))
    #     except fit.InversionError:
    #         logging.warn("Covariance matrix failed, skipping this tmin,tmax {},{}".format(*tmin))
    #     else:
    #         print "fit using ranges {}".format(tmin)
    #         break
    # print "done"
    # exit()

    # ranges = [(10,25), (10,25)]
    multicor = mergecors(cors, zip(*ranges)[0], zip(*ranges)[1])
    funct = functions[options.function](Nt=options.period, ranges=ranges)
    fit.fit(funct, multicor, min(multicor.times), max(multicor.times),
            filestub=options.output_stub, return_chi=False, return_quality=True, options=options)
    print "done"


    logging.info("starting fit with mergedcorrelator")


def plot(cors, funct, time_starts, time_ends, averages, stds, chi, options):
    import plot_helpers
    import newton
    dt = 1

    fig = plt.figure()

    corplots = {}
    corplots[0] = plt.subplot(221)
    corplots[1] = plt.subplot(222)
    emassplots = {}
    emassplots[0] = plt.subplot(223)
    emassplots[1] = plt.subplot(224)


    for i in [0,1]:
        corplots[i].errorbar(cors[i].times,
                             cors[i].average_sub_vev().values(),
                             yerr=cors[i].jackknifed_errors().values(), fmt='o')

        emass = cors[i].periodic_effective_mass(dt, fast=False, period=options.period)
        emass_errs = cors[i].periodic_effective_mass_errors(dt, fast=False).values()
        emassplots[i].errorbar(emass.keys(), emass.values(),
                               yerr=emass_errs, fmt='o')

    fn = funct.individual(Nt=options.period)
    for i in [0, 1]:
        #X = np.array(cors[i].times)
        X = np.linspace(time_starts[i], time_ends[i])
        X = np.arange(time_starts[i], time_ends[i]+1, 1)
        Xe = np.arange(time_starts[i], time_ends[i], 1)

        Y = fn.formula((averages[0], averages[i+1]), X)
        Ydt = fn.formula((averages[0], averages[i+1]), X+dt)
        corplots[i].plot(X, Y)
        guess = (1.0 / float(dt)) * np.log(Y/Ydt)
        fitemass = {}
        Yd = dict(enumerate(Y, start=time_starts[i]))
        guess = dict(enumerate(guess, start=time_starts[i]))
        for t in range(time_starts[i], time_ends[i]):
            fitemass[t] = newton.newton_cosh_for_m(t, t+dt, Yd, guess[t], options.period)
        emassplots[i].plot(fitemass.keys(),fitemass.values())

    if options.output_stub:
        filestub = options.output_stub
        logging.info("Saving plot to {}".format(filestub+".png"))
        fig.set_size_inches(18.5, 10.5)
        plt.savefig(filestub+".png")
    else:
        plt.show()

if __name__ == "__main__":

    # add the fit parser, but then override to use many correlators.
    # confict handler is required for that
    parser = argparse.ArgumentParser(description="compute fits", parents=[fitparser], conflict_handler='resolve')


    parser.add_argument("-i", "--inputfile", action='append', type=str, required=True,
                           help="Correlator files to read from")
    parser.add_argument("-ts", "--time-start", type=int, action='append', required=False,
                            help="first time slice to start a fit, can be a list of times")
    parser.add_argument("-te", "--time-end", type=int, action='append', required=False,
                            help="last time slice to fit, can be a list of times")

    parser.add_argument("--debug_heavycorrection", action="store_true",
                        help="make the correlator symmetric")


    function_list = inspect.getmembers(sys.modules["simul_fitfunctions"], inspect.isclass)
    functions = {name: f for name, f in function_list}
    parser.add_argument("-f", "--function", choices=functions.keys(),
                           required=True, default="twocor_periodic_exp", help="function to fit to")

    args = parser.parse_args()


    if args.verbose:
        logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.DEBUG)
        logging.error("Wtf")
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


    if args.output_stub and args.skip_done:
        filename = args.output_stub+".boot"
        if args.jackknife:
            filename = args.output_stub+".jack"
        try:
            if os.stat(filename).st_size > 0:
                logging.info(".boot file exists and not empty, skip fit")
                logging.info("{}".format(filename))
                exit(0)
            else:
                logging.warn(".boot file exists but is empty!")
        except OSError:
            logging.info("running fit")

    consistant_argument_counts(args)

    cors = []
    for i in range(len(args.inputfile)):
        corrfile = args.inputfile[i]

        vev1 = args.vev[i]
        vev2 = vev1
        if args.vev2:
            vev2 = args.vev2[i]

        cor = build_corr.corr_and_vev_from_pickle(corrfile, vev1, vev2)


        corsym = cor.determine_symmetry()
        if corsym is None:
            logging.error("called with symmetric but correlator isnt")
            raise RuntimeError("called with symmetric but correlator isnt")
        logging.info("correlator found to be {}".format(corsym))
        cor.make_symmetric()
        cor.prune_invalid(delete=True, sigma=args.prune)

        if args.debug_heavycorrection:
            import heavy_correction_dicts
            beta = [b for b in ["4.17", "4.35", "4.47"] if b in args.inputfile[i]][0]
            mtype = [b for b in ["m0", "m1", "m2", "m3", "m4", "m5"] if b in args.inputfile[i]][0]
            print beta
            print mtype
            d = heavy_correction_dicts.get_heavy_correction(beta, mtype)
            cor.multiply_by_value_dict(d)


        if args.bin:
            cor = cor.reduce_to_bins(args.bin)


        cors.append(cor)

    logging.info("fitting using start times {} and end times{} ".format( args.time_start, args.time_end))

    DONE = False

    if not args.time_start:
        auto_fit(cors, options=args)
    else:
        while not DONE:
            multicor = mergecors(cors, args.time_start, args.time_end)
            multicor.symmetric = True
            multicor.symmetry = "symmetric" # hack
            multicor.period = args.period

            funct = functions[args.function](Nt=args.period, ranges=zip(args.time_start, args.time_end))

            funct.stride = args.tstride

            logging.info("starting fit with mergedcorrelator")
            try:
                averages, stds, chi = fit.fit(funct, multicor,
                                              min(multicor.times), max(multicor.times), bootstraps=args.bootstraps,
                                              filestub=args.output_stub, return_chi=True,
                                              return_quality=False, writecor=False, tstride=args.tstride, options=args)
            except (fit.InversionError, InvalidFit):
                logging.error("Could not invert, trying larger stride time {}->{}".format(args.tstride, args.tstride+1))
                args.tstride += 1
                funct.stride = args.tstride
                for s,e in zip(args.time_start, args.time_end):
                    if (e-s)/args.tstride < 4:
                        raise RuntimeError("Shrunk too far still cant invert")
                continue

            DONE = True

            if args.plot:
                plot(cors, funct, args.time_start, args.time_end, averages, stds, chi, args)
