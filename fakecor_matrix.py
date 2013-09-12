#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
import logging
import correlator
import build_corr
import pylab
import argparse
import os
import zfactor
import fakecor

from fitfunctions import *  # noqa
import inspect
import sys




if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Generate a fake correlator")
    parser.add_argument("-o", "--output_stub", type=str, required=True,
                        help="stub of name to write output to")
    parser.add_argument("-cfg", "--configs", type=int, required=False, default=50,
                        help="number of configs to make")
    parser.add_argument("-t", "--times", type=int, required=False, default=20,
                        help="number of times")
    parser.add_argument("-z", "--overlapsfile", type=str, required=True,
                        help="file to read the zfactors")
    parser.add_argument("-m", "--masses", type=float, nargs="+", required=True,
                        help="list of masses for each level")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="increase output verbosity")

    args = parser.parse_args()
    if args.verbose:
        logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.DEBUG)
        logging.debug("Verbose debuging mode activated")
    else:
        logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.INFO)

    N = len(args.masses)

    raw_Zs = zfactor.read_coeffs_file(args.overlapsfile)
    assert len(raw_Zs)==N*N, "Length missmatch in levels and coeffs"
    print raw_Zs
    Zs = np.matrix(raw_Zs.identities.values.reshape((N, N))).T
    print Zs
    print "blah"
    for i in range(N):
        for j in range(N):
            print "Correlator_{}{}".format(i,j)
            amps = [Zs[i,level]*Zs[j,level] for level in range(N)]
            print "newaps", amps
            cor = fakecor.make_fake_cor(args.configs, args.times, amps, args.masses)
            cor.writefullfile(args.output_stub+"{}_{}".format(i,j), comp=True)

    logging.info("done")
