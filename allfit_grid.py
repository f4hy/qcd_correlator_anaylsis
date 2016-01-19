#!/usr/bin/env python2
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.lines as mlines
import matplotlib.patches as patches
import logging                  # Including many defaults, can be removed if unneeded
import argparse
import os
import numpy as np
import pandas as pd
import math
import re

def allfit_grid(options):
    """ Plot a grid fr different fit ranges """
    logging.debug("Called with {}".format(options))

    data = {}
    data["chi"] = np.full((32,32), np.inf)
    data["chidof"] = np.full((32,32), np.inf)


    tmin = 100
    tmax = -1

    with open(options.inputfile) as readfile:
        for line in readfile:
            if "t=" in line:

                sline = line.split("=")[-1].split()
                ts = int(sline[0])
                te = int(sline[-1].strip("\n-"))
                tmin = min(tmin,ts)
                tmax = max(tmax,te)
                readfile.next()
                for dataline in readfile:
                    if dataline.startswith("-"):
                        chiline = readfile.next()
                        schiline = chiline.split(",")
                        chi = float(schiline[0].split("=")[-1])
                        chidof = float(schiline[1].split("=")[-1])
                        qual = float(schiline[2].split(" ")[-1])
                        data["chi"][ts][te] = chi
                        data["chidof"][ts][te] = chidof
                        break
                    else:
                        sdata = dataline.split()
                        name = sdata[0]
                        average = float(sdata[2])
                        std = float(sdata[4])
                        if name not in data.keys():
                            data[name] = np.full((32,32), np.inf)
                            data[name+"_std"] = np.full((32,32), np.inf)
                            data[name+"_rel"] = np.full((32,32), np.inf)
                        data[name][ts][te] = average
                        data[name+"_std"][ts][te] = std
                        data[name+"_rel"][ts][te] = std/average

    fig, axe = plt.subplots(1)

    cmap = mpl.cm.cool
    cmap = mpl.cm.spring
    cmap = mpl.cm.rainbow
    #cmap = mpl.colors.ListedColormap(['blue','black','red'])
    #cmap = mpl.cm.cool
    bounds=[-1,0,0.5,1.0]
    #norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

    if options.data in data.keys():
        plotdata = data[options.data]
        if options.data == "chidof":
            tile = np.percentile(plotdata[plotdata<1], 50)
            plotdata[plotdata > 1.0] = np.inf
    else:
        logging.error("do not have data on {} options are {}".format(options.data, data.keys()))

    img = plt.imshow(plotdata, interpolation='nearest', cmap=cmap)

    plt.xlim(tmin-0.5,tmax+0.5)
    plt.ylim(tmin-0.5,tmax+0.5)

    # make a color bar
    plt.colorbar(img,cmap=cmap)

    plt.show()


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Plot a grid showing data for differnet fit ragnes")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="increase output verbosity")
    parser.add_argument("-o", "--output_stub", type=str, required=False,
                        help="stub of name to write output to")
    parser.add_argument("-i", "--inputfile", type=str, required=True,
                        help="name of file to write read")
    parser.add_argument("-d", "--data", type=str, required=True,
                        help="what data to plot")
    parser.add_argument('--err', nargs='?', type=argparse.FileType('w'),
                        default=None)
    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.DEBUG)
        logging.debug("Verbose debuging mode activated")
    else:
        logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.INFO)

    if args.err is not None:
        root = logging.getLogger()
        ch = logging.StreamHandler(args.err)
        ch.setLevel(logging.ERROR)
        formatter = logging.Formatter('%(levelname)s: %(message)s')
        ch.setFormatter(formatter)
        root.addHandler(ch)


    allfit_grid(args)
