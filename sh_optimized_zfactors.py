#!/usr/bin/env python2
import argparse
import logging
import plot_files
import numpy as np
import pandas as pd
import pandas_reader
import re

from level_identifier import readops
from level_identifier import read_file
from zfactor import read_coeffs_file


class overlaps(object):
    "class for rotation coeffs and zfactors"
    def __init__(self, panda_data):
        self.data = panda_data
        self.ops, self.levels = map(int,re.search("(\d+)0+(\d+)", str(self.data.index[-1])).groups(1))
        logging.info("created object with {}ops and {}levels".format(self.ops,self.levels))

    def index_format(self, op, level):
        if 0 <= level <= self.levels and 0 <= op <= self.ops:
            return int("{}{:03d}".format(op,level))
        else:
            raise IndexError("Level or op out of bounds")

    def get_entry(self, op, level):
        i = self.index_format(op,level)
        print i
        return self.get_index(i)

    def get_index(self, i):
        return self.data.ix[i].identities

    def get_level(self, level):
        print [self.get_entry(op,level) for op in range(1,self.ops+1)]

    def get_op(self, op):
        print [self.get_entry(op,level) for level in range(1,self.levels+1)]


        # print self.data


def sh_optimized_zfacts():
    shops = readops(args.single_hadron_ops)
    mhops = readops(args.full_hadron_ops)
    fullz = overlaps(read_file(args.full_zfactors))
    rotco = overlaps(read_coeffs_file(args.rotation_coeffs))
    print shops
    indicies = [mhops.index(o) for o in shops]
    print indicies
    # print fullz
    # print rotco
    print [mhops[i] for i in indicies]

    print fullz.index_format(1,1)
    print fullz.index_format(103,82)

    print fullz.get_entry(1,1)
    # print fullz.get_index(1001)
    print fullz.get_level(60)
    print fullz.get_op(64)

    # for i in range(1,105):
    #     index = int("{}001".format(i))
    #     print index
    #     print fullz.get_index(index)



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="compute the SH optimized overlaps from two diags")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="increase output verbosity")
    parser.add_argument("-o", "--output_stub", type=str, required=False,
                        help="stub of name to write output to")
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
