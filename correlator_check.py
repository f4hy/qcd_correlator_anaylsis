#!/usr/bin/env python2
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import logging
import argparse
import numpy as np
import build_corr

from scipy import stats

def check(f, options):

    cor = build_corr.corr_and_vev_from_pickle(f, None, None)

    logging.info("correlator with {} cfgs and {} times ".format(len(cor.configs),len(cor.times)))
    logging.info("correlator with period {}".format(cor.period))
    logging.info("correlator with period {}".format(cor.period_check(96)))

    if options.symmetry:
        corsym = cor.determine_symmetry(recheck=True)
        if "PP" in f:
            if corsym != "symmetric":
                logging.error("PP correlator not symmetric")
                exit(-1)
        if "A4P" in f or "PA4" in f:
            if corsym != "anti-symmetric":
                logging.error("A4P correlator not anti-symmetric: {}".format(f))
                exit(-1)
    # sanity_check(data)

    if options.rel_err:
        c = cor.average_sub_vev()
        e = cor.jackknifed_errors()
        #print c
        logging.info("max jk rel error {}".format(max(e[t] / c[t] for t in cor.times)))

    if options.stats:
        e = cor.jackknifed_errors()
        for t in cor.times:
            ct = np.array(cor.get(time=t).values())
            print ct.mean(), np.median(ct)
            print ct.std(), e[t]
            upperquart = stats.scoreatpercentile(ct, 75)
            lowerquart = stats.scoreatpercentile(ct, 25)
            print upperquart, ct.mean()+ct.std()
            print lowerquart, ct.mean()-ct.std()

            exit(-1)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="check a correlator file")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="increase output verbosity")
    parser.add_argument("--stats", action="store_true",
                        help="check statistics")
    parser.add_argument("-s", "--symmetry", action="store_true",
                        help="check symmetry")
    parser.add_argument("-r", "--rel_err", action="store_true",
                        help="check relative error")
    parser.add_argument("-o", "--output-stub", type=str, required=False,
                        help="stub of name to write output to")
    parser.add_argument('files', metavar='f', type=str, nargs='+',
                        help='files to plot')
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


    for f in args.files:

        check(f, args)
