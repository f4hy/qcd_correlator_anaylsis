#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
import logging
import argparse
import plot_files
import pandas as pd
from operator_tranlator import translate
import math

OUTPUT = 25
logging.addLevelName(OUTPUT, "OUTPUT")

def readops(filename):
    with open(filename) as opfile:
        raw_operators = opfile.readlines()
        if len(raw_operators) == 1:
            operators = raw_operators[0].split()
        else:
            operators = map(str.strip, raw_operators)
    return operators


def read_file(filename):
    txt = plot_files.lines_without_comments(filename)
    filetype = plot_files.determine_type(txt)
    if filetype == "paren_complex":
        parse_pair = plot_files.parse_pair
        df = pd.read_csv(txt, delimiter=' ', names=["id", "identities", "error"],
                         converters={1: parse_pair, 2: parse_pair}, index_col=0)
    else:
        df = pd.read_csv(txt, delimiter=' ', names=["id", "identities", "error"], index_col=0)
    return df


def make_bar_plot(inputfile, cols, output_stub, mode, ns, opnames=None, maxplots=1000):
    def square(x):
        return np.absolute(x)**2
    df = read_file(inputfile).apply(square)
    ops = list(set([i/1000 for i in df.index]))
    levels = list(set([(i % 1000) for i in df.index]))
    N = int(np.sqrt(len(df)))
    largest_zfactor = max(df.identities)
    untranslated = opnames
    if opnames:
        opnames = [translate(l) for l in opnames]
    plots = ops
    if mode == "level":
        plots = levels

    plot = {}
    # plt.figure(figsize=(10, 6))
    plt.tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='off')
    plt.rcParams['xtick.major.size'] = 0

    rows = min(int(math.ceil(float(len(plots))/cols)), int(math.ceil(float(maxplots)/cols)))
    if not args.seperate:
        f, layout = plt.subplots(nrows=rows, ncols=cols)
    # f.set_size_inches(19.2, 12.0)
    # f.set_dpi(100)
    for plot_index in plots[:maxplots]:
        if args.seperate:
            f, ax = plt.subplots(1)
        i = (plot_index-1)/cols
        j = (plot_index-1) % cols
        logging.info("making plot {} {} {}".format(plot_index, i, j))
        if not args.seperate:
            if len(plots[:maxplots]) <= cols or cols == 1:
                ax=layout[j]
            else:
                ax=layout[i][j]

        indexes = [plot_index+op*1000 for op in ops]
        if j > 0 and not args.seperate:
            ax.set_yticks([])
        else:
            ax.set_ylabel("$|Z|^2$", fontweight='bold')
            # ax.set_yticks([0, 1])

        if mode == "ops":
            with open(args.ordering) as orderfile:
                ordering = [int(i.strip())+1 for i in orderfile.readlines()] # orderings have levels indexed from 0, but these are indexed from 1

            indexes = [plot_index*1000+level for level in ordering]
            ticklabelpad = plt.rcParams['xtick.major.pad']
            ax.set_xticks(np.array(levels)+1.5)
            ax.set_xticklabels([n if n % 5 == 0 else "" for n in levels])
            ax.annotate('Level', xy=(1, 0), xytext=(-10, -ticklabelpad*2.2), ha='left', va='top',
                        xycoords='axes fraction', textcoords='offset points')
            if opnames:
                ax.set_title("{}".format(opnames[plot_index-1]))
            else:
                ax.set_title("operator {}".format(plot_index))
        else:
            ax.set_title("level {}".format(plot_index-1))
            ticklabelpad = plt.rcParams['xtick.major.pad']
            ax.annotate('Operator', xy=(1, 0), xytext=(-10, -ticklabelpad*2.2), ha='left', va='top',
                        xycoords='axes fraction', textcoords='offset points')

        values = df.ix[indexes].identities.values
        errors = df.ix[indexes].error.values
        if args.no_errors:
            errors = None
        # plots[(i, j)] = ax.bar(ops, values, yerr=errors)
        if mode == "level":
            plot[(i, j)] = ax.bar(ops[:ns], values[:ns], 1, color="g", yerr=errors)
            plot[(i, j)] = ax.bar(ops[ns:], values[ns:], 1, color="b", yerr=errors)
            ax.set_xticks(np.array(ops)+0.5)
            ax.set_xticklabels([o if o % 5 == 0 else "" for o in ops])
            plt.title("level {}".format(plot_index))
        else:
            color = "b"
            if plot_index <= ns:
                color = "g"
            plot[(i, j)] = ax.bar(levels, values, 1, color=color, yerr=errors, ecolor="r")
        ax.set_ylim([0, np.ceil(largest_zfactor)])
        ax.set_xlim(xmin=1)
        if args.seperate:
            ax.set_ylim([0, max(values)+max(errors)])
            outfilename = output_stub+"{}.png".format(untranslated[plot_index-1])
            logging.info("Saving plot to {}".format(outfilename))
            plt.savefig(outfilename)
            plt.clf()

    # end plot loop

    if args.seperate:
        return

    if args.title:
        f.suptitle(args.title)

    if(output_stub):
        plt.rcParams.update({'font.size': 8})
        plt.tight_layout(pad=2.0, h_pad=1.0, w_pad=2.0)
        logging.info("Saving plot to {}".format(output_stub+".png"))
        plt.savefig(output_stub+".png")
        # logging.info("Saving plot to {}".format(output_stub+".eps"))
        # plt.savefig(output_stub+".eps")
    else:
        plt.tight_layout()
        plt.show()


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Plot the coffeficents of diagonalized operators")
    parser.add_argument("-i", "--inputfile", type=str, required=True,
                        help="Correlator file to read from")
    parser.add_argument("-o", "--output_stub", type=str, required=False,
                        help="stub of name to write output to")
    parser.add_argument("--ordering", type=str, required=True,
                        help="file which contains the ordering")
    parser.add_argument("-c", "--columns", type=int, required=False,
                        help="number of columns to make the plot", default=2)
    parser.add_argument("-m", "--mode", type=str, required=False, choices=["level", "ops"],
                        help="select the mode of how to plot the levels", default="ops")
    parser.add_argument("-ns", "--number-singlehadrons", type=int, required=False, default=0,
                        help="Number of single hadrons to distinguish")
    parser.add_argument("-mx", "--max-plots", type=int, required=False, default=100,
                        help="Maximum number of plots")
    parser.add_argument("-n", "--names", type=str, required=False,
                        help="operator names file")
    parser.add_argument("-t", "--title", type=str, required=False,
                        help="plot title", default=None)
    parser.add_argument("-ne", "--no-errors", action="store_true", default=None,
                        help="suppress the error bars")
    parser.add_argument("-sep", "--seperate", action="store_true", default=None,
                        help="put each on their own plot")
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
        names = readops(args.names)
        print names

    make_bar_plot(args.inputfile, args.columns, args.output_stub, args.mode, args.number_singlehadrons, names, maxplots=args.max_plots)
