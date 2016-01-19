#!/usr/bin/env python2
import argparse
import logging
import plot_files
import numpy as np
import pandas as pd
import pandas_reader
import build_corr

def writeformat(x):
    """Format complex number into paren grouped format"""
    return [int(x[0]) ,"({},{})".format(np.real(x[1]), np.imag(x[1]))]

def all_same(items):
    return all(x == items[0] for x in items)

def ratio_files(files, options):
    # data = [plot_files.read_file(f) for f in files]
    c1 = build_corr.corr_and_vev_from_files(files[0], None, None)
    c2 = build_corr.corr_and_vev_from_files(files[1], None, None)

    a1 = c1.average_sub_vev()
    a2 = c2.average_sub_vev()

    ratio = {}
    for t in c1.times:
        print t
        print a1[t], a2[t]
        print a2[t]/a1[t]
        ratio[t] = a2[t]/a1[t]

    if options.output_stub:
        filename = options.output_stub+".ratio"
        logging.info("writing ratio to {}".format(filename))
        with open(filename, 'w') as ofile:
            for t,v in ratio.iteritems():
                ofile.write("{}, {}\n".format(t,v))


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

    ratio_files(args.files, args)
