#!/usr/bin/env python
import matplotlib.pyplot as plt
import logging
import argparse
from plot_files import read_file
import determine_operators


def plot_matrix(input_directory, fileformat, operators):
    N = len(operators)
    plots = {}
    for i, row in enumerate(operators):
        for j, col in enumerate(operators):
            print fileformat.format(row, col)
            filename = fileformat.format(row, col)
            ax = plt.subplot2grid((N, N), (i, j))
            df = read_file(input_directory+filename)
            plots[(i, j)] = ax.errorbar(df.time.values, df.correlator.values, yerr=df.error.values,
                                        linestyle="none", marker="o", label="{}_{}".format(row, col))
            # ax.legend(fancybox=True, shadow=True)
            if i == 0:
                ax.set_title(col)
            if j == 0:
                ax.set_ylabel(row)

    plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot matrix of files")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="increase output verbosity")
    parser.add_argument("-f", "--format", type=str, required=True,
                        help="fromat of the correlator files in the directory\n\n"
                        "e.g. {}-{}.A1gp.conn.dat where {} are replaced with operator strings")
    parser.add_argument("-r", "--operators", action='append', required=False,
                        help="operator to make \n\n e.g. -r a1pp_0_optype0_op1")
    parser.add_argument("-i", "--input-dir", type=str, required=True,
                        help="directory to read files from")
    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.DEBUG)
        logging.debug("Verbose debuging mode activated")
    else:
        logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.INFO)

    if not args.operators:
        logging.info("Operators not specified, attempting to automagically determine")
        ops = determine_operators.matching_operators(args.input_dir, args.format)
        print ops
        if not ops:
            print "Error: no operators found"
            parser.print_help()
            parser.exit()
        args.operators = ops

    plot_matrix(args.input_dir, args.format, args.operators)
