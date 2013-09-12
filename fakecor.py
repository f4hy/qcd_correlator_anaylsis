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





def make_fake_cor(cfgs, times, amps, masses):
    assert(len(amps) == len(masses))
    cfgs = list(range(cfgs))
    times = list(range(times))
    data = {}
    vev = {}
    for c in cfgs:
        vev[c] = 0.0
        tmp = {}
        for t in times:
            tmp[t] = sum([amps[i] * np.exp((-1.0*masses[i]) * t) for i in range(len(amps))])
            tmp[t] += np.random.normal(0,0.01,1) * tmp[t]
            data[c] = tmp

    return correlator.Correlator.fromDataDicts(data, vev, vev)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Generate a fake correlator")
    parser.add_argument("-o", "--output_stub", type=str, required=True,
                        help="stub of name to write output to")
    parser.add_argument("-cfg", "--configs", type=int, required=False, default=50,
                        help="number of configs to make")
    parser.add_argument("-t", "--times", type=int, required=False, default=20,
                        help="number of times")
    parser.add_argument("-a", "--amps", type=float, nargs="+", required=True,
                        help="amplitudes for each exp")
    parser.add_argument("-m", "--masses", type=float, nargs="+", required=True,
                        help="masses for each exp")
    parser.add_argument("-b", "--amp2", type=float, required=False, default=0.0,
                        help="amplitude for second exp")
    parser.add_argument("-m2", "--mass2", type=float, required=False, default=1.5,
                        help="mass for second exp")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="increase output verbosity")

    args = parser.parse_args()
    if args.verbose:
        logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.DEBUG)
        logging.debug("Verbose debuging mode activated")
    else:
        logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.INFO)


    print args.amps
    print args.masses
    cor = make_fake_cor(args.configs, args.times, args.amps, args.masses)
    cor.writefullfile(args.output_stub, comp=True)
