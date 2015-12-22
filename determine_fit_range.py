#!/usr/bin/env python2
import logging                  # Including many defaults, can be removed if unneeded
import argparse
import os
import numpy as np
import pandas as pd
import math
import re
import glob
import build_corr


def determine_fit_range(options):
    """ get fit range """
    logging.debug("Called with {}".format(options))

    files = glob.glob("SymDW_sHtTanh_b2.0_smr3_*x{}x*b{}*0/data/*{}*{}*{}-{}*{}".format(options.period, options.beta, options.heavy, options.flavor, options.smearing1, options.smearing2, options.operator))

    if len(files) < 1:
        logging.error("no files")
        exit(-1)

    print files
    files = files
    #fhandles = [open(f) for f in files]
    maxt = (options.period/2)-3
    rel_err = {t:0 for t in range(1,options.period/2)}

    cors = [build_corr.corr_and_vev_from_files(f, None, None) for f in files]

    for c in cors:
        if options.operator == "PP":
            c.make_symmetric()
        if options.operator == "A4P":
            c.make_symmetric(anti=True)
            maxt = (options.period/2)-8

    for cor in cors:
        asv = cor.average_sub_vev()
        errs = cor.jackknifed_errors()
        masserrs = cor.effective_mass_errors(1)
        for t in range(4,options.period/2):
            rel_err[t] = max(errs[t]/asv[t], rel_err[t])
            if t > options.period/2 and masserrs[t] > masserrs[t-2]*2:
                logging.warning("Emass error doubled in two slices!")
                print t, masserrs[t] , masserrs[t-2]
                maxt = min(maxt, t-2)


    print rel_err
    for t,e in rel_err.iteritems():
        if t > maxt:
            break
        if e > 0.15:
            print e
            maxt = t-4
            break

    maxt = maxt
    mint = -1

    for cor in cors:
        emass = cor.cosh_effective_mass(1, fast=True, period=options.period)
        errs = cor.cosh_effective_mass_errors(1, fast=True, period=options.period)
        # emass = cor.effective_mass(1)
        # errs = cor.effective_mass_errors(1)
        prev = emass[3]
        for t in range(4, maxt):
            if emass[t] > prev:
                print "increased at", t
                mint = max(mint,t)
                break
            if emass[t]+errs[t] > emass[t-1] and emass[t-1]+errs[t-1] > emass[t-2] and emass[t]+errs[t] > emass[t-2]:
                print "const at", t
                mint = max(mint,t)
                break

            prev = emass[t]

    if mint < 0:
        logging.error("could not find const")
        exit(-1)

    params = [options.period, options.beta, options.heavy, options.flavor, options.smearing1, options.smearing2, options.operator]
    params = map(str,params)
    string = "_".join(params)
    ofilename = options.output_stub+"/"+string
    logging.info("fit range deteremined to be {} {}".format(mint,maxt))
    logging.info("writing fit range to {}".format(ofilename))

    with open(ofilename, 'w') as outfile:
        outfile.write("{} {}\n".format(mint, maxt))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="~/develop/anaysis_iroiro/")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="increase output verbosity")
    parser.add_argument("-o", "--output_stub", type=str, required=False,
                        help="stub of name to write output to")
    parser.add_argument('period', type=int, help='PERIOD')
    parser.add_argument('operator', type=str, help='OPERATOR')
    parser.add_argument('heavy', type=str, help='heavy')
    parser.add_argument('flavor', type=str, help='FLAVOR')
    parser.add_argument('smearing1', type=int, help='smearing1')
    parser.add_argument('smearing2', type=int, help='smearing2')
    parser.add_argument('beta', type=str, help='beta')
    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.DEBUG)
        logging.debug("Verbose debuging mode activated")
    else:
        logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.INFO)

    determine_fit_range(args)
