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
    print raw_Zs
    Zs = np.matrix(raw_Zs.identities.values.reshape((N, N))).T
    print Zs
    print "blah"
    for i in range(N):
        for j in range(N):
            print "Correlator_{}{}".format(i,j)
            Ai = Zs[i,0]
            Bi = Zs[i,1]
            Aj = Zs[j,0]
            Bj = Zs[j,1]
            amp1 = Ai*np.conj(Aj)
            amp2 = Bi*np.conj(Bj)
            print Ai, Aj, Bi, Bj
            # exit()
            cor = fakecor.make_fake_cor(args.configs, args.times, amp1, args.masses[0],
                                        amp2, args.masses[1])
            cor.writefullfile(args.output_stub+"{}_{}".format(i,j), comp=True)

    logging.info("done")
