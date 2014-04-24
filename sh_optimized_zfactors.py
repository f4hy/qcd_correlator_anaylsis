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



    N = len(OptZ.keys())
    Ncols=args.columns
    rows = int(math.ceil(float(N)/Ncols))
    fig, ax = plt.subplots(ncols=Ncols, nrows=rows)
    # print OptZ.keys(), range(1,N+1)
    for m in range(1,N+1):
        i = (m-1)/Ncols
        j = (m-1) % Ncols
        ax[i][j].bar(range(len(OptZ[m])), OptZ[m], 1.0, color="b")
        ax[i][j].set_title("SH-opt level{}".format(m))
        ax[i][j].set_ylim((0,max(OptZ[m])))
        ax[i][j].set_ylabel("$|Z|^2$", fontweight='bold')
        ax[i][j].set_xlabel("Level", fontweight='bold')


    # plt.ylim((0,np.max(OptZ.values())))


    if args.output_stub:
        logging.info("Saving shopt_zfactors to {}".format(args.output_stub+".out"))
        with open(args.output_stub+".out", 'w') as outfile:
            for level,d in OptZ.iteritems():
                outfile.write("{}, {}\n".format(level, ", ".join(map(str,d))))
        logging.info("Saving shcandidates to {}".format(args.output_stub+".singleresonances"))
        with open(args.output_stub+".singleresonances", 'w') as resfile:
            for level,d in OptZ.iteritems():
                resfile.write("{}, {}\n".format(level, np.argmax(d)))

        plt.rcParams.update({'font.size': 10})
        fig.set_size_inches(18.5,10.5)
        plt.tight_layout()
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
    parser.add_argument("-c", "--columns", type=int, default=3, required=False,
                        help="Number of columns")
    parser.add_argument("-z", "--full_zfactors", type=str, required=True,
                        help="zfactors")
    parser.add_argument("-s", "--single_hadron_ops", type=str, required=True,
                        help="list of single hadron ops")
    parser.add_argument("-m", "--full_hadron_ops", type=str, required=True,
                        help="full list of ops")
    parser.add_argument("-r", "--rotation_coeffs", type=str, required=True,
                        help="single hadtron rotation coeffs")
    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.DEBUG)
        logging.debug("Verbose debuging mode activated")
    else:
        logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.INFO)

    sh_optimized_zfacts()
