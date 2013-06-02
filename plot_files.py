#!/usr/bin/env python
import matplotlib.pyplot as plt
from matplotlib import mpl
import pandas_reader
import logging
import build_corr
import argparse
import fit
from matplotlib.widgets import CheckButtons
from compiler.ast import flatten
import os
import pandas as pd

from fitfunctions import *  # noqa
import inspect
import sys
import StringIO
import re

pair = re.compile(r'\(([^,\)]+),([^,\)]+)\)')


def parse_pair(s):
    print "testing"
    print s
    if s:
        return complex(*map(float, pair.match(s).groups()))
    else:
        return ""


def read_file(f):
    # txt = f.read()
    # df = pd.read_csv(f, delimiter=' ',  names=["correlator", "error"], index_col=0, converters={1: parse_pair, 2: parse_pair})

    # for line in f:
    #     print line
    #     if line[0] == "#":
    #         continue
    #     print line.split(" ")
    #     t,c,e = line.split(" ")
    #     print c
    #     print parse_pair(c)
    #     print e
    #     print parse_pair(e)
    #     exit()
    # df = pd.read_csv(f, delimiter=' ', comment="#", names=["time", "correlator", "error"],
    #                  converters={1: parse_pair, 2: parse_pair})

    df = pd.read_csv(f, delimiter=' ', comment="#", names=["time", "correlator", "error"])
    # df = pd.read_csv(f, delimiter=',', comment="#", names=["time", "correlator", "error"])

    print df.head(20)
    print np.real(df.correlator[1])
    print parse_pair(df.correlator[1])
    exit()
    return df


def plot_files(files):
    markers = ['o', "D" , "^", "<", ">", "v", "x", "p", "8"]
    # colors, white sucks
    colors = [c for c in mpl.colors.colorConverter.colors.keys() if c!='w']
    print files
    plots = {}
    for index,filename in enumerate(files):
        f = open(filename, "r")
        basename = os.path.basename(filename)
        df = read_file(f)
        print df.head(20)
        print df.time.values, df.correlator.values, df.error.values
        plots[basename] = plt.errorbar(df.time.values+(index*0.1), df.correlator.values, yerr=df.error.values, linestyle="none", c=colors[index%len(colors)], marker=markers[index%len(markers)])
    print len(colors), len(markers)


    def toggle_errorbar_vis(ebarplot):
        for i in flatten(ebarplot):
            if i:
                i.set_visible(not i.get_visible())


    def func(label):
        toggle_errorbar_vis(plots[label])
        plt.draw()


    rax = plt.axes([0.85, 0.8, 0.1, 0.15])
    check = CheckButtons(rax, [os.path.basename(f) for f in args.files], [True]*len(plots))
    check.on_clicked(func)

    plt.show()



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="compute fits")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="increase output verbosity")
    # parser.add_argument('files', metavar='f', type=argparse.FileType('r'), nargs='+',
    #                     help='files to plot')
    parser.add_argument('files', metavar='f', type=str, nargs='+',
                        help='files to plot')
    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.DEBUG)
        logging.debug("Verbose debuging mode activated")
    else:
        logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.INFO)

    plot_files(args.files)
