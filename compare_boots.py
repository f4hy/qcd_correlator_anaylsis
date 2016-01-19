#!/usr/bin/env python2
import logging
import argparse
import os
import numpy as np
import pandas as pd
import math
import re


def compare_boots(options):
    """ Compare two data files to find the differences """
    logging.debug("Called with {}".format(options))

    col_names = ["bootstrap", "mass", "amp1", "amp2", "amp3", "amp4"]


    df1 = pd.read_csv(options.files[0], comment='#', names=col_names, index_col=0)
    df2 = pd.read_csv(options.files[1], comment='#', names=col_names, index_col=0)


    df1=df1.dropna(axis=1,how='all')
    df2=df2.dropna(axis=1,how='all')


    (df1.mean() - df2.mean()) > min(df1.std())



    stats = df1.mean() - df2.mean()
    stats = pd.DataFrame(stats, columns=["difference"])
    stats["percent diff"] = 100*abs(df1.mean() - df2.mean()) / (0.5*(df1.mean() + df2.mean()))
    stats["sigma1"] = abs(df1.mean() - df2.mean()) / df1.std()
    stats["sigma2"] = abs(df1.mean() - df2.mean()) / df2.std()
    stats["error percent diff"] = 100*abs(df1.std() - df2.std()) / (0.5*(df1.std() + df2.std()))
    # stats["sigma"] = abs(df1.mean() - df2.mean()) / df2.std()

    logging.info("\n{}".format(stats))

    if any(stats["sigma1"] > 1.0):
        logging.error("not within one sigma")
    else:
        logging.info("within one sigma")

    if all(df1.std() < df2.std()):
        logging.info("file2 has larger errors")
        if options.force_improvement:
            logging.error("file2 is not better than file1")
    elif all(df1.std() > df2.std()):
        logging.info("file1 has larger errors")
    else:
        logging.warn("mixed error sizes")
        logging.warn(df1.std())
        logging.warn(df2.std())



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compare two bootstraped fits")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="increase output verbosity")
    parser.add_argument("--force_improvement", action="store_true",
                        help="increase output verbosity")
    parser.add_argument("-o", "--output_stub", type=str, required=False,
                        help="stub of name to write output to")
    parser.add_argument('files', metavar='f', type=str, nargs='+',
                        help='input files')
    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.DEBUG)
        logging.debug("Verbose debuging mode activated")
    else:
        logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.INFO)

    if args.output_stub is not None:
        root = logging.getLogger()
        errfilename = args.output_stub+".err"
        errfilehandler = logging.FileHandler(errfilename, delay=True)
        errfilehandler.setLevel(logging.ERROR)
        formatter = logging.Formatter('%(levelname)s: %(message)s')
        errfilehandler.setFormatter(formatter)
        root.addHandler(errfilehandler)
        logfilename = args.output_stub+".log"
        logfilehandler = logging.FileHandler(logfilename, delay=True)
        logfilehandler.setLevel(logging.INFO)
        formatter = logging.Formatter('%(levelname)s: %(message)s')
        logfilehandler.setFormatter(formatter)
        root.addHandler(logfilehandler)


    compare_boots(args)
