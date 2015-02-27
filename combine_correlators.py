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

def writerealformat(x):
    """Format complex number into paren grouped format"""
    return [int(x[0]) , float(x[1])]

def all_same(items):
    return all(x == items[0] for x in items)


def read_full_correlator(filename, emass=False):
    try:
        cor = build_corr.corr_and_vev_from_files_pandas(filename, None, None)
    except AttributeError:
        logging.info("Failed to read with pandas, reading normal")
        cor = build_corr.corr_and_vev_from_files(filename, None, None)

    if emass:
        emasses = cor.cosh_effective_mass(1)
        times = emasses.keys()
        data = [emasses[t] for t in times]
        errors = [cor.cosh_effective_mass_errors(1)[t]  for t in times]
    else:
        times = cor.times
        data = [cor.average_sub_vev()[t] for t in times]
        errors = [cor.jackknifed_errors()[t]  for t in times]

    d = {"time": times, "correlator": data, "error": errors, "quality": [float('NaN') for t in times]}
    df = pd.DataFrame(d)
    return df


def combine_files(files, options):
    data = [plot_files.read_file(f) for f in files]
    assert all_same([df.shape for df in data]), "They are not all the same size!"
    if options.function == "average":
        new = sum(data)/len(data)
    elif options.function == "ratio":
        if len(data) != 2:
            raise RuntimeError("Ratio requires exactly 2 correlators")
        new = data[0][["time", "correlator"]]
        new["correlator"] = new["correlator"] / data[1]["correlator"]
    elif options.function == "sum":
        new = sum(data)


    formated = new[['time','correlator']]
    formated["time"] = formated["time"].astype(int)
    print formated
    formated.columns = ["#time", "correlator"]
    print formated
    outputtext = formated.to_csv(None, sep=",", index=False, header=False, index_label="#time").replace(",", ", ")
    print outputtext
    if options.output_stub:
        outfile = "{}.dat".format(options.output_stub)
        with open(outfile, 'w') as ofile:
            ofile.write(outputtext)
            logging.info("wrote ave cor to {}".format(outfile))
    else:
        print outputtext

if __name__ == "__main__":
    functs = ["average", "ratio", "sum"]
    parser = argparse.ArgumentParser(description="average data files")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="increase output verbosity")
    parser.add_argument("-o", "--output_stub", type=str, required=False,
                        help="stub of name to write output to")
    parser.add_argument("-f", "--function", type=str, required=False, choices=functs,
                        help="function to use to combine the data")
    parser.add_argument('files', metavar='f', type=str, nargs='+',
                        help='files to plot')
    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.DEBUG)
        logging.debug("Verbose debuging mode activated")
    else:
        logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.INFO)

    combine_files(args.files, args)