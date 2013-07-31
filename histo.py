#!/usr/bin/env python
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from matplotlib import mpl
import logging
import argparse
from matplotlib.widgets import CheckButtons
from compiler.ast import flatten
import os
import pandas as pd
import pandas_reader
from operator_tranlator import translate
import math

from fitfunctions import *  # noqa
from cStringIO import StringIO

import re

pair = re.compile(r'\(([^,\)]+),([^,\)]+)\)')


def make_histogram(data, numbins=10, norm=False):

    print "norm", norm

    fig = plt.figure()
    if np.iscomplexobj(data):
        realplot = plt.subplot(211)
        realplot.set_title("real")
        realdata = np.real(data)
        _, bins, _ = realplot.hist(realdata, numbins, normed=norm)
        if norm:
            bincenters = 0.5*(bins[1:]+bins[:-1])
            y = mlab.normpdf( bincenters, np.mean(realdata), np.std(realdata))
            realplot.plot(bincenters, y, 'r--', linewidth=1)


        imagplot = plt.subplot(212)
        imagplot.set_title("imag")
        imagdata = np.imag(data)
        _, bins, _ = imagplot.hist(np.imag(data), numbins, facecolor="green", normed=norm)
        if norm:
            bincenters = 0.5*(bins[1:]+bins[:-1])
            y = mlab.normpdf( bincenters, np.mean(imagdata), np.std(imagdata))
            imagplot.plot(bincenters, y, 'r--', linewidth=1)
    else:
        plt.hist(data, bins)
        if norm:
            bincenters = 0.5*(bins[1:]+bins[:-1])
            y = mlab.normpdf( bincenters, np.mean(data), np.std(data))
            plt.plot(bincenters, y, 'r--', linewidth=1)

    plt.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="plot a histogram of a file for a single time")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="increase output verbosity")
    parser.add_argument("-t", "--time", type=int, required=True,
                        help="time slice to histogram", default=None)
    parser.add_argument("-b", "--bins", type=int, required=False, default=10,
                        help="number of bins for the histogram")
    parser.add_argument("-n", "--norm", action="store_true",
                        help="normalize and draw normal distribution")
    parser.add_argument('datafile', metavar='f', type=str, help='file to plot')
    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.DEBUG)
        logging.debug("Verbose debuging mode activated")
    else:
        logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.INFO)

    data = pandas_reader.read_single_time_paraenformat(args.datafile, args.time)


    make_histogram(data, args.bins, args.norm)
