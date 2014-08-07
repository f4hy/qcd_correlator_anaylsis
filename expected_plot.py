#!/usr/bin/env python2
import matplotlib.pyplot as plt
import itertools
import logging
import argparse
import re
import numpy as np
from primary_operators import read_expected_levels
import irreps

def main():
    if args.verbose:
        logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.DEBUG)
        logging.debug("Verbose debuging mode activated")
    else:
        logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.INFO)

    logging.info("Reading expected levels from {}".format(args.inputfile))
    levels = read_expected_levels("T1up_1")
    plotlevels(levels)

def read_expected_levels(channel):
    logging.info("opening expected levels from {}".format(args.inputfile))
    expected_level_file = open(args.inputfile, "r")
    chan = channel.split("_")[0]
    start = False
    expectedleveltxt = ""
    for line in expected_level_file:
        if " " + chan in line:
            start = True
        if start:
            if "(" in line:
                expectedleveltxt += line
            if not line.strip():
                break

    levellines = [level.strip().split() for level in expectedleveltxt.split("\n") if level]
    levelpairs = [ (float(i[0]), particle_type(i[-1])) for i in levellines ]
    return levelpairs

def particle_type(s):
    """ takes an expected level and returns a generic particle type"""
    imap = {0: "isoscalar" , "1/2": "K", 1: "pi"}
    underscores = s.count("_")
    if underscores is 0:
        return "single"
    if underscores is 5:
        ss = s.split("_")
        return "-".join(imap[irreps.particle_name(ss[i]).I] for i in (0,2))
    if underscores is 7:
        return "3hadron"
    if underscores is 9:
        return "4hadron"
    raise NotImplementedError("Unsupported level, only supports mesons for now")

def plotlevels(levels):
    colors = ['b', 'r', 'k', 'm', 'c', 'y', 'g']
    markers = ['o', "D", "^", "<", ">", "v"]
    cmarkers = zip(colors*5, markers*5)
    cmarkers.reverse()
    colordefs = {"single": "".join(cmarkers.pop())}
    values = [i[0] for i in levels]
    types = [i[1] for i in levels]
    for t in types:
        if t not in colordefs.keys():
            colordefs[t] = "".join(cmarkers.pop())

    fig, ax = plt.subplots()
    plots = []
    for t in colordefs.keys():
        indexes = [i for i,s in enumerate(types) if s==t]
        y = [values[i] for i in indexes]
        plots.append(ax.plot(indexes, y, colordefs[t], label=t, ms=12))

    leg = plt.legend(fancybox=True, shadow=True, loc=2)
    plt.xlim(-1,len(values)+2)
    plt.ylim(0,3)
    ax.set_xlabel("Expected Level", fontweight='bold', fontsize=30)
    ax.set_ylabel("m", fontweight='bold', fontsize=30)

    if(args.output_stub):
        fig.set_size_inches(18.5, 10.5)
        plt.rcParams.update({'font.size': 20})
        logging.info("Saving plot to {}".format(args.output_stub+".png"))
        plt.savefig(args.output_stub+".png", dpi=100)
        # logging.info("Saving plot to {}".format(output_stub+".eps"))
        # plt.savefig(output_stub+".eps")
        return
    plt.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="plot a set of data files")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="increase output verbosity")
    parser.add_argument("-o", "--output-stub", type=str, required=False,
                        help="stub of name to write output to")
    parser.add_argument("-y", "--yrange", type=float, required=False, nargs=2,
                        help="set the yrange of the plot", default=None)
    parser.add_argument("-3", "--threshold", type=float, required=False,
                        help="Draw a line where 3 particle threshold is")
    parser.add_argument("-i", "--inputfile", type=str, required=True,
                        help="expected level file to read from")
    args = parser.parse_args()

    main()
