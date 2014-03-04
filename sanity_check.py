#!/usr/bin/env python
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import logging
import argparse
import pandas_reader
import pandas
import numpy as np

from scipy import stats

def sanity_check(data):

    data = data.apply(np.real,reduce=False)
    statistics = pandas.DataFrame()
    statistics["mean"] = data.mean(1)
    statistics["median"] = data.median(1)
    statistics["std"] = data.std(1)
    statistics["skew"] = statistics["mean"] - statistics["median"]
    statistics["percentskew"] = (statistics["mean"] - statistics["median"])/statistics["median"]
    statistics.to_csv(args.output_stub+".stats",sep="\t", float_format="%10.5g")

    biggest_skew = max(statistics["percentskew"])
    logging.info("largest skew: {}".format(biggest_skew))
    if biggest_skew > 50.0:
        logging.error("Data in {} is highly skewed!!".format(args.datafile))

    outliers = pandas.DataFrame()
    outliers_big = {}
    outliers_small = {}
    outlierfile = open(args.output_stub+".outliers", "w")
    outliercount = {}
    outlierfile.write( "Time | outliers too large  | outliers too small\n")
    for time in data.index:
        upperquart = stats.scoreatpercentile(data.loc[time], 75)
        lowerquart = stats.scoreatpercentile(data.loc[time], 25)
        iqr = upperquart-lowerquart
        sigma=5
        upperbound = float(np.real(statistics["median"][time]+iqr*sigma))
        lowerbound = float(np.real(statistics["median"][time]-iqr*sigma))
        outlierfile.write("{:<5}".format(time))
        toobig = data.loc[time][(data.loc[time] > upperbound)].index.values
        if toobig.size > 0:
            outliers_big[time] = "".join(toobig)
            outlierfile.write("{:<30}".format("".join(toobig).replace("correlator","")))
            for o in [s.replace("correlator","config") for s in toobig]:
                if o in outliercount:
                    outliercount[o] += 1
                else:
                    outliercount[o] = 1
        else:
            outlierfile.write( "                              ",)
        toosmall = data.loc[time][(data.loc[time] < lowerbound)].index.values
        if toosmall.size > 0:
            outliers_small[time] = "".join(toosmall)
            outlierfile.write(" {:<30}\n".format("".join(toosmall).replace("correlator","")))
            for o in [s.replace("correlator","config") for s in toosmall]:
                if o in outliercount:
                    outliercount[o] += 1
                else:
                    outliercount[o] = 1
        else:
            outlierfile.write("\n")
    for config,count in outliercount.items():
        outlierfile.write("{} has outliers on {} of {} times\n".format(config,count,len(data.index)))
        if count == len(data.index):
            logging.error("Config {} has outliers on all times!!".format(config))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="plot a histogram of a file for a single time")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="increase output verbosity")
    parser.add_argument("-o", "--output-stub", type=str, required=True,
                        help="stub of name to write output to")
    parser.add_argument('datafile', metavar='f', type=str, help='file to plot')
    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.DEBUG)
        logging.debug("Verbose debuging mode activated")
    else:
        logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.INFO)

    data = pandas_reader.read_configcols_paraenformat(args.datafile)

    sanity_check(data)
