#!/usr/bin/env python
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import logging
import argparse
import pandas_reader
import plot_files

from fitfunctions import *  # noqa

import re

pair = re.compile(r'\(([^,\)]+),([^,\)]+)\)')


def make_histogram(data, options, output_stub, numbins=10, norm=False):

    print "norm", norm

    if options.logarithm:
        data = np.log(data)

    fig = plt.figure()          # noqa
    if np.iscomplexobj(data):
        realplot = plt.subplot(211)
        realplot.set_title("real")
        realdata = np.real(data)
        _, bins, _ = realplot.hist(realdata, numbins, normed=norm)
        if norm:
            bincenters = 0.5*(bins[1:]+bins[:-1])
            y = mlab.normpdf(bincenters, np.mean(realdata), np.std(realdata))
            realplot.plot(bincenters, y, 'r--', linewidth=1)

        imagplot = plt.subplot(212)
        imagplot.set_title("imag")
        imagdata = np.imag(data)
        _, bins, _ = imagplot.hist(np.imag(data), numbins, facecolor="green", normed=norm)
        if norm:
            bincenters = 0.5*(bins[1:]+bins[:-1])
            y = mlab.normpdf(bincenters, np.mean(imagdata), np.std(imagdata))
            imagplot.plot(bincenters, y, 'r--', linewidth=1)
    else:
        _, bins, _ = plt.hist(data, numbins, normed=norm)
        if norm:
            bincenters = 0.5*(bins[1:]+bins[:-1])
            y = mlab.normpdf(bincenters, np.mean(data), np.std(data))
            plt.plot(bincenters, y, 'r--', linewidth=1)
            plt.yticks([])

    # if options.title:
    #     plt.suptitle(options.title)

    if options.zero:
        plt.xlim(0, plt.xlim()[1])

    if(output_stub):
        logging.info("Saving plot to {}".format(output_stub+".png"))
        plt.savefig(output_stub+".png",dpi=200)
    else:
        plt.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="plot a histogram of a file for a single time")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="increase output verbosity")
    parser.add_argument("-l", "--logarithm", action="store_true",
                        help="take the log of the data first")
    parser.add_argument("-t", "--time", type=int, required=False,
                        help="time slice to histogram", default=None)
    parser.add_argument("-b", "--bins", type=int, required=False, default=100,
                        help="number of bins for the histogram")
    parser.add_argument("-n", "--norm", action="store_true",
                        help="normalize and draw normal distribution")
    parser.add_argument("-c", "--column", type=int, required=False,
                        help="column to histogram", default=None)
    parser.add_argument("-o", "--output-stub", type=str, required=False,
                        help="stub of name to write output to")
    parser.add_argument("--title", type=str, required=False,
                        help="title to add to plot")
    parser.add_argument("-z", "--zero", action="store_true",
                        help="include zero")
    parser.add_argument('datafile', metavar='f', type=str, help='file to plot')
    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.DEBUG)
        logging.debug("Verbose debuging mode activated")
    else:
        logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.INFO)

    if args.time:
        #data = pandas_reader.read_single_time_paraenformat(args.datafile, args.time)
        data = pandas_reader.read_single_time_commaformat(args.datafile, args.time)
    elif args.column is not None:
        data =plot_files.read_file(args.datafile, list(range(100)))[[args.column]].values
        if all(np.isnan(data)):
            logging.error("Column is all nan")
            logging.shutdown()
            exit(-1)
    else:
        with open(args.datafile) as dataf:
            txt = [line.strip(",") for line in dataf.read().split() if not line.startswith("#")]
            data= map(float,txt)

    logging.info("making histogram of data of {} points".format(len(data)))
    make_histogram(data, args, args.output_stub, args.bins, args.norm)
