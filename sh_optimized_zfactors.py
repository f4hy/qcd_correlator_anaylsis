#!/usr/bin/env python2
import argparse
import logging
import plot_files
import numpy as np
import pandas as pd
import pandas_reader
import re
import matplotlib.pyplot as plt
import math

from level_identifier import readops
from level_identifier import read_file
from zfactor import read_coeffs_file
from plot_helpers import ncols


class overlaps(object):
    "class for rotation coeffs and zfactors"
    def __init__(self, panda_data):
        self.data = panda_data
        self.ops, self.levels = map(int,re.search("(\d+)0+(\d\d+)", str(self.data.index[-1])).groups(1))
        logging.info("created object with {}ops and {}levels".format(self.ops,self.levels))

    def index_format(self, op, level):
        if 0 <= level <= self.levels and 0 <= op <= self.ops:
            return int("{}{:03d}".format(op,level))
        else:
            raise IndexError("Level or op out of bounds")

    def get_entry(self, op, level):
        i = self.index_format(op,level)
        return self.get_index(i)

    def get_index(self, i):
        return self.data.ix[i].identities

    def get_level(self, level):
        return [self.get_entry(op,level) for op in range(1,self.ops+1)]

    def get_op(self, op):
        return [self.get_entry(op,level) for level in range(1,self.levels+1)]


        # print self.data


def sh_optimized_zfacts():
    shops = readops(args.single_hadron_ops)
    mhops = readops(args.full_hadron_ops)
    rotco = overlaps(read_coeffs_file(args.rotation_coeffs))
    fullz = overlaps(read_file(args.full_zfactors))
    indicies = [mhops.index(o)+1 for o in shops]

    OptZ = {}
    # Main calculation below
    # TODO this should be done as matrix multiplication, would be WAY faster but this works.
    for m in range(1,rotco.levels+1):
        value = np.array([np.abs(np.array(np.matrix( np.conj(rotco.get_level(m))) * np.matrix([fullz.get_entry(i,l) for i in indicies]).T))**2 for l in range(1,fullz.levels+1)]).flatten()
        if np.all(value==0):
            break
        OptZ[m] = value

    with open(args.ordering) as orderfile:
        ordering = [int(i.strip()) for i in orderfile.readlines()]


    N = len(OptZ.keys())
    if args.number:
        N = min(N,args.number)
    if args.columns:
        Ncols = args.columns
    else:
        Ncols = ncols(N)
    rows = int(math.ceil(float(N)/Ncols))
    if not args.seperate:
        fig, ax = plt.subplots(ncols=Ncols, nrows=rows)
    with open(args.SHordering) as SHorderfile:
        SHordering = [int(i.strip()) for i in SHorderfile.readlines()]
    for index, m in enumerate(SHordering):
        if index >= N:
            break
        i = (index)/Ncols
        j = (index) % Ncols
        if args.seperate:
            fig, axe = plt.subplots(1)
        else:
            if N <= Ncols:
                axe=ax[j]
            else:
                axe=ax[i][j]
        reordered = [OptZ[m+1][reorderedlevel] for reorderedlevel in ordering]

        axe.bar(range(len(reordered)), reordered, 1.0, color="b")
        # axe.set_title("SH-opt level{}".format(index))
        axe.set_ylim((0,max(reordered)))
        axe.set_ylabel("$|Z|^2$", fontweight='bold', fontsize=30)
        axe.set_xlabel("Level", fontweight='bold', fontsize=30)
        axe.tick_params(axis='both', which='major', labelsize=20)
        plt.tight_layout()
        if args.seperate:
            outfilename = args.output_stub+"_{}.eps".format(index)
            logging.info("Saving plot to {}".format(outfilename))
            plt.savefig(outfilename)
            plt.clf()

    if args.seperate:
        return

    if args.output_stub:
        logging.info("Saving unreordered shopt_zfactors to {}".format(args.output_stub+".nonreordered.out"))
        with open(args.output_stub+".out", 'w') as outfile:
            for level,d in OptZ.iteritems():
                outfile.write("{}, {}\n".format(level, ", ".join(map(str,d))))
        logging.info("Saving shcandidates to {}".format(args.output_stub+".singleresonances"))
        with open(args.output_stub+".singleresonances", 'w') as resfile:
            for level,d in OptZ.iteritems():
                fixd = np.nan_to_num(d)
                if max(fixd) > 0:
                    resfile.write("{}, {}\n".format(level, np.argmax(fixd)))

        plt.rcParams.update({'font.size': 10})
        fig.set_size_inches(18.5,8.5)
        plt.tight_layout()
        if args.eps:
            logging.info("Saving plot to {}".format(args.output_stub+".eps"))
            plt.savefig(args.output_stub+".eps")
        else:
            logging.info("Saving plot to {}".format(args.output_stub+".png"))
            plt.savefig(args.output_stub+".png",dpi=100)
    else:
        plt.tight_layout()
        plt.show()




if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="compute the SH optimized overlaps from two diags")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="increase output verbosity")
    parser.add_argument("-o", "--output_stub", type=str, required=False,
                        help="stub of name to write output to")
    parser.add_argument("--ordering", type=str, required=True,
                        help="file which contains the ordering")
    parser.add_argument("--SHordering", type=str, required=True,
                        help="file which contains the ordering")
    parser.add_argument("-c", "--columns", type=int, required=False,
                        help="Number of columns")
    parser.add_argument("-z", "--full_zfactors", type=str, required=True,
                        help="zfactors")
    parser.add_argument("-s", "--single_hadron_ops", type=str, required=True,
                        help="list of single hadron ops")
    parser.add_argument("-m", "--full_hadron_ops", type=str, required=True,
                        help="full list of ops")
    parser.add_argument("-r", "--rotation_coeffs", type=str, required=True,
                        help="single hadtron rotation coeffs")
    parser.add_argument("-N", "--number", type=int, required=False,
                        help="max number of plots")
    parser.add_argument("-sep", "--seperate", action="store_true", default=None,
                        help="put each on their own plot")
    parser.add_argument("--eps", action="store_true",
                        help="plot eps instead of png")
    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.DEBUG)
        logging.debug("Verbose debuging mode activated")
    else:
        logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.INFO)

    sh_optimized_zfacts()
