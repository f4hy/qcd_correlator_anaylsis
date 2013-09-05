#!/usr/bin/env python
import numpy as np
import logging
import argparse
import plot_files
import pandas as pd
import re

pair = re.compile(r'\(([^,\)]+),([^,\)]+)\)')


def parse_pair(s):
    if s:
        return complex(*map(float, pair.match(s).groups()))
    else:
        return ""


def read_file(filename):
    txt = plot_files.lines_without_comments(filename)
    df = pd.read_csv(txt, delimiter=' ', names=["time", "correlator", "error"],
                     converters={1: parse_pair, 2: parse_pair}, skipinitialspace=True, index_col=0)
    return df


def compared_cormatrix(corwild, reconwild, ops):
    chi_sqr = 0.0
    for col, src in enumerate(ops):
        for row, snk in enumerate(ops):
            logging.debug("Reading snk:{}, src:{}".format(snk, src))
            raw_c = read_file(corwild.format(snk, src))
            raw_rc = read_file(reconwild.format(snk, src))
            diff = raw_c.correlator-raw_rc.correlator
            mag = np.absolute(diff)
            sqr = np.square(mag)
            summed = np.sum(sqr)
            chi_sqr += summed
    logging.info("chi_sqr is {}".format(chi_sqr))


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Compute the Zfactors from the correlator and diagonalized coeffs")
    parser.add_argument("-ic", "--inputcorrelatorformat", type=str, required=True,
                        help="Correlator wild card to read from")
    parser.add_argument("-ir", "--inputreconstructed", type=str, required=True,
                        help="Reconstructed correlator wild card to read from")
    parser.add_argument("-ops", "--operators", type=str, nargs="+", required=True,
                        help="operator strings, order matters!")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="increase output verbosity")

    args = parser.parse_args()
    if args.verbose:
        logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.DEBUG)
        logging.debug("Verbose debuging mode activated")
    else:
        logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.INFO)

    compared_cormatrix(args.inputcorrelatorformat, args.inputreconstructed, args.operators)
