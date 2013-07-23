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



def make_bar_plot(inputfile, cols, output_stub, mode, ns, opnames=None):
    df = read_file(inputfile).apply(np.absolute)

    ops = list(set([i/1000 for i in df.index]))
    levels = list(set([i % 1000 for i in df.index]))

    N = int(np.sqrt(len(df)))
    largest_zfactor = max(df.identities)

    opnames = ["".join(n.split("-")[1:]) if n.startswith("iso") else n for n in opnames ]
    plots = ops
    if mode == "level":
        plots = levels

    print "levels", levels
    print "ops", ops
    print "plots", plots
    plot = {}
    plt.figure(figsize=(10, 6))
    plt.suptitle("Zfactors by {}".format(mode))
    for plot_index in plots:
        print "making plot {}".format(plot_index)
        i = (plot_index-1)/cols
        j = (plot_index-1) % cols

        ax = plt.subplot2grid((len(plots)/cols+1, cols), (i, j))
        indexes = [plot_index+op*1000 for op in ops]
        if j > 0:
            ax.set_yticks([])
        else:
            plt.ylabel("Z", fontweight='bold', rotation='horizontal')
            ax.set_yticks([0,1])

        plt.xlabel("operator")
        if mode == "ops":
            indexes = [plot_index*1000+level for level in levels]
            plt.xlabel("                  level", fontweight='bold')
            if opnames:
                plt.title("{}".format(opnames[plot_index-1]))
            else:
                plt.title("operator {}".format(plot_index))
        values = df.ix[indexes].identities.values
        # errors = df.ix[indexes].error.values
        # plots[(i, j)] = ax.bar(ops, values, yerr=errors)
        if mode == "level":
            plot[(i, j)] = ax.bar(ops[:ns], values[:ns], 1, color="g")
            plot[(i, j)] = ax.bar(ops[ns:], values[ns:], 1, color="b")
            ax.set_xticks(np.array(ops)+0.5)
            ax.set_xticklabels([o if o % 5 == 0 else "" for o in ops])
            plt.title("level {}".format(plot_index))
        else:
            color = "b"
            if plot_index <= ns:
                color = "g"
            plot[(i, j)] = ax.bar(levels, values, 1, color=color, lw=2)
            ax.set_xticks(np.array(levels)+0.5)
            ax.set_xticklabels([n if n % 5 == 0 else "" for n in levels])
        plt.ylim([0, np.ceil(largest_zfactor)])
        plt.xlim(xmin=1)
    # end plot loop

    plt.tight_layout()
    plt.subplots_adjust(wspace=0.2, hspace=1.0)
    if(output_stub):
        plt.rcParams.update({'font.size': 12})
        logging.info("Saving plot to {}".format(output_stub+".png"))
        plt.savefig(output_stub+".png", dpi=300)
    else:
        plt.rcParams.update({'font.size': 14})
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
    parser.add_argument("-n", "--names", type=str, required=False,
                        help="operator names file")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="increase output verbosity")

    args = parser.parse_args()
    if args.verbose:
        logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.DEBUG)
        logging.debug("Verbose debuging mode activated")
    else:
        logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.INFO)

    names = None
    if args.names:
        with open(args.names) as namefile:
            names = namefile.readline().split()
        print names

    make_bar_plot(args.inputfile, args.columns, args.output_stub, args.mode, args.number_singlehadrons, names)
