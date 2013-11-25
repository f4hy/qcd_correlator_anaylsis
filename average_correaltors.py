#!/usr/bin/env python2
import argparse
import logging
import plot_files
import numpy as np
import pandas as pd
import pandas_reader

def writeformat(x):
    """Format complex number into paren grouped format"""
    return [int(x[0]) ,"({},{})".format(np.real(x[1]), np.imag(x[1]))]


def average_files(files, outfile):
    data = [plot_files.read_file(f) for f in files]
    r = pandas_reader.read_normal_paraenformat
    data = [r(f) for f in files]
    ave = sum(data)/len(data)
    formated = ave.apply(writeformat,axis=1)
    formated.columns = ["#time", "correlator"]
    print formated
    formated.to_csv(outfile, sep=" ", index=False, index_label="#time")
    logging.info("wrote ave cor to {}".format(outfile))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="average data files")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="increase output verbosity")
    parser.add_argument("-o", "--output_stub", type=str, required=False,
                        help="stub of name to write output to")
    parser.add_argument('files', metavar='f', type=str, nargs='+',
                        help='files to plot')
    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.DEBUG)
        logging.debug("Verbose debuging mode activated")
    else:
        logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.INFO)

    average_files(args.files, args.output_stub)
