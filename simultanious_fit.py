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
            cor = build_corr.corr_and_vev_from_files_pandas(corrfile, vev1, vev2)
            cors.append(cor)
        except AttributeError:
            logging.info("Failed to read with pandas, reading normal")
            cor = build_corr.corr_and_vev_from_files(corrfile, vev1, vev2)
        if args.symmetric:
            cor.make_symmetric()
        if args.antisymmetric:
            cor.make_symmetric(anti=True)

        cors.append(cor)


    multicor = mergecors(cors, args.time_start, args.time_end)
    funct = functions[args.function](Nt=args.period, ranges=zip(args.time_start, args.time_end))

    logging.info("starting fit with mergedcorrelator")
    fit.fit(funct, multicor, min(multicor.times), max(multicor.times), options=args)
