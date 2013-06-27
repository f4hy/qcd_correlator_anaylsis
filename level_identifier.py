#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
import logging
import argparse
import plot_files
import pandas as pd

OUTPUT = 25
logging.addLevelName(OUTPUT, "OUTPUT")


def read_file(filename):
    txt = plot_files.lines_without_comments(filename)
    # parse_pair = plot_files.parse_pair
    # df = pd.read_csv(txt, delimiter=' ', names=["id", "identities", "error"],
    #                  converters={1: parse_pair, 2: parse_pair}, index_col=0)
    df = pd.read_csv(txt, delimiter=' ', names=["id", "identities", "error"], index_col=0)
    return df



def make_bar_plot(inputfile, cols, output_stub, mode, ns):
    df = read_file(inputfile).apply(np.absolute)

    N = int(np.sqrt(len(df)))
    ops = range(1, N+1)

    plots = {}
    plt.figure(figsize=(10, 6))
    for op1 in ops:
        i = (op1-1)/cols
        j = (op1-1) % cols

        ax = plt.subplot2grid((N/cols+1, cols), (i, j))
        indexes = [op1+op2*1000 for op2 in ops]
        plt.xlabel("operator")
        if mode == "ops":
            indexes = [op1*1000+op2 for op2 in ops]
            plt.xlabel("levels")
        values = df.ix[indexes].identities.values
        # errors = df.ix[indexes].error.values
        # plots[(i, j)] = ax.bar(ops, values, yerr=errors)
        if mode == "level":
            plots[(i, j)] = ax.bar(ops[:ns], values[:ns], color="g")
            plots[(i, j)] = ax.bar(ops[ns:], values[ns:], color="b")
        else:
            color = "b"
            if op1 <= ns:
                color = "g"
            plots[(i, j)] = ax.bar(ops, values, color=color)


        plt.ylim([0, 1])

    plt.tight_layout()
    if(output_stub):
        plt.rcParams.update({'font.size': 5})
        logging.info("Saving plot to {}".format(output_stub+".png"))
        plt.savefig(output_stub+".png", dpi=200)
    else:
        plt.show()


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Plot the coffeficents of diagonalized operators")
    parser.add_argument("-i", "--inputfile", type=str, required=True,
                        help="Correlator file to read from")
    parser.add_argument("-o", "--output_stub", type=str, required=False,
                        help="stub of name to write output to")
    parser.add_argument("-c", "--columns", type=int, required=False,
                        help="number of columns to make the plot", default=2)
    parser.add_argument("-m", "--mode", type=str, required=False, choices=["level", "ops"],
                        help="select the mode of how to plot the levels", default="level")
    parser.add_argument("-ns", "--number-singlehadrons", type=int, required=False, default=0,
                        help="Number of single hadrons to distinguish")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="increase output verbosity")

    args = parser.parse_args()
    if args.verbose:
        logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.DEBUG)
        logging.debug("Verbose debuging mode activated")
    else:
        logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.INFO)

    make_bar_plot(args.inputfile, args.columns, args.output_stub, args.mode, args.number_singlehadrons)
