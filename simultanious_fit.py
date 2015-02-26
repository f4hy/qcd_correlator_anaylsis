#!/usr/bin/env python
import fit
import inspect
import sys
from simul_fitfunctions import *  # noqa
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
    for conf in configs:
        dummyvev[conf] = 0.0
        dataindex = 0
        confdata = {}
        for i in range(0,len(cors)):
            for t in range(tmins[i], tmaxs[i]+1):
                confdata[dataindex] = cors[i].get(config=conf, time=t)
                dataindex +=1
        newdata[conf] = confdata


    return correlator.Correlator.fromDataDicts(newdata, dummyvev, dummyvev)


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

    consistant_argument_counts(args)

    cors = []
    for i in range(len(args.inputfile)):
        corrfile = args.inputfile[i]

        vev1 = args.vev[i]
        vev2 = vev1
        if args.vev2:
            vev2 = args.vev2[i]

        try:
            logging.info("reading {}".format(corrfile))
            cor = build_corr.corr_and_vev_from_files_pandas(corrfile, vev1, vev2)
            cors.append(cor)
        except AttributeError:
            logging.info("Failed to read with pandas, reading normal")
            cor = build_corr.corr_and_vev_from_files(corrfile, vev1, vev2)

        print corrfile
        if "A4P" in corrfile:
            cor.make_symmetric(anti=True)
        elif "PP" in corrfile:
            cor.make_symmetric()


        cors.append(cor)


    if not args.time_start:
        auto_fit(cors, options=args)
    else:
        multicor = mergecors(cors, args.time_start, args.time_end)
        funct = functions[args.function](Nt=args.period, ranges=zip(args.time_start, args.time_end))

        logging.info("starting fit with mergedcorrelator")
        fit.fit(funct, multicor, min(multicor.times), max(multicor.times), filestub=args.output_stub, return_chi=False, return_quality=True, options=args)
