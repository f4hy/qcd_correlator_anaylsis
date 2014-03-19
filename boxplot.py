#!/usr/bin/env python2
import matplotlib.pyplot as plt
import matplotlib as mpl
import logging
import argparse
from matplotlib.widgets import CheckButtons
from compiler.ast import flatten
import os
import pandas as pd
from operator_tranlator import translate
import math

from fitfunctions import *  # noqa
from cStringIO import StringIO

import re

pair = re.compile(r'\(([^,\)]+),([^,\)]+)\)')


def parse_pair(s):
    if s:
        return complex(*map(float, pair.match(s).groups()))
    else:
        return ""


def lines_without_comments(filename, comment="#"):
    s = StringIO()
    with open(filename) as f:
        for line in f:
            if not line.startswith(comment):
                s.write(line)
    s.seek(0)
    return s

def myconverter(s):
    try:
        return np.float(s.strip(','))
    except:
        return np.nan

def removecomma(s):
    return int(s.strip(','))


def determine_type(txt):
    firstline = txt.readline()
    txt.seek(0)
    if "(" in firstline and ")" in firstline:
        return "paren_complex"
    if "," in firstline:
        return "comma"
    return "space_seperated"


def read_file(filename):
    with open(filename, 'r') as f:
        first_line = f.readline()
    names=[s.strip(" #") for s in first_line.split(",")[0:-2]]
    print names
    txt = lines_without_comments(filename)
    filetype = determine_type(txt)
    if filetype == "paren_complex":
        df = pd.read_csv(txt, delimiter=' ', names=names,
                         converters={1: parse_pair, 2: parse_pair})
    if filetype == "comma":
        df = pd.read_csv(txt, sep=",", delimiter=",", names=names, skipinitialspace=True, delim_whitespace=True, converters={0: removecomma, 1: myconverter})
    if filetype == "space_seperated":
        df = pd.read_csv(txt, delimiter=' ', names=names)
    return df


def get_fit(filename):
    with open(filename) as f:
        for line in f:
            if "fit" in line:
                logging.info("found fit info: {}".format(line))
                fitrange = re.search("\(([0-9]+),([0-9]+)\)", line)
                tmin, tmax = int(fitrange.group(1)), int(fitrange.group(2))
                mass = float(re.search("m=(.*?) ", line).group(1))
                error = float(re.search("e=(.*?) ", line).group(1))
                qual = re.search("qual:(.*)", line).group(1)
                return (tmin, tmax, mass, error, qual)
    raise RuntimeError("No fit info")


def label_names_from_filelist(filelist):
    names = filelist
    basenames = [os.path.basename(filename) for filename in names]
    names = basenames
    if any(basenames.count(x) > 1 for x in basenames):
        logging.debug("two files with same basename, cant use basenames")
        names = filelist

    if len(names) < 2:
        return names
    prefix = os.path.commonprefix(names)
    if len(prefix) > 1:
        names = [n[len(prefix):] for n in names]
    postfix = os.path.commonprefix([n[::-1] for n in names])
    if len(postfix) > 1:
        names = [n[:-len(postfix)] for n in names]
    names = [(n.strip(" _") if n != "" else "base") for n in names]
    return names


def add_fit_info(filename, ax=None):
    if not ax:
        ax = plt

    try:
        tmin, tmax, mass, error, quality = get_fit(filename)
        ax.plot(range(tmin, tmax+1), [mass]*len(range(tmin, tmax+1)))
        ax.plot(range(tmin, tmax+1), [mass+error]*len(range(tmin, tmax+1)), ls="dashed", color="b")
        ax.plot(range(tmin, tmax+1), [mass-error]*len(range(tmin, tmax+1)), ls="dashed", color="b")
        digits = -1.0*round(math.log10(error))
        formated_error = int(round(error * (10**(digits + 1))))
        formated_mass = "{m:.{d}}".format(d=int(digits) + 1, m=mass)
        # ax.annotate("{m}({e})".format(m=formated_mass, e=formated_error), xy=(tmax,mass),
        #              xytext=(tmax+1, mass+error))
        return "{m}({e}) qual:{q:.4}".format(m=formated_mass, e=formated_error, q=quality)
    except RuntimeError:
        logging.error("File {} had no fit into".format(filename))


def boxplot_files():
    markers = ['o', "D", "^", "<", ">", "v", "x", "p", "8"]
    # colors, white sucks
    colors = [c for c in mpl.colors.colorConverter.colors.keys() if c != 'w' and c != "g"]
    plots = {}
    tmin_plot = {}
    has_colorbar = False
    labels = label_names_from_filelist(args.files)
    #labels = [translate(l) for l in labels]
    layout = None
    data = []
    f, ax = plt.subplots()
    for index, label, filename in zip(range(len(args.files)), labels, args.files):

        # for index, label in enumerate(labels):
        mark = markers[index % len(markers)]
        color = colors[index % len(colors)]
        df = read_file(filename)

        if args.seperate:
            data.append(df.mass.values)
        else:
            plots[label] = plt.boxplot(df.mass.values, widths=0.5, patch_artist=True)
            plt.setp(plots[label]["boxes"], color=color)
            plt.setp(plots[label]["whiskers"], color=color)
            plt.setp(plots[label]["fliers"], color=color)

    if args.seperate:
        plt.boxplot(data, widths=0.5, patch_artist=True)
        xticknames = plt.setp(ax, xticklabels=labels)
        plt.setp(xticknames, rotation=45, fontsize=8)
    if not args.seperate:
        leg = plt.legend(fancybox=True, shadow=True)
    # else:
    #     plt.tight_layout(pad=0.0, h_pad=0.0, w_pad=0.0)
    #     if has_colorbar:
    #         f.subplots_adjust(right=0.95)
    #         cbar_ax = f.add_axes([0.96, 0.05, 0.01, 0.9])
    #         f.colorbar(tmin_plot[label], cax=cbar_ax)

    if(args.output_stub):
        if args.title:
            f.suptitle(args.title)
        f.set_size_inches(18.5,10.5)
        plt.rcParams.update({'font.size': 5})
        #plt.tight_layout(pad=2.0, h_pad=1.0, w_pad=2.0)
        plt.tight_layout()
        logging.info("Saving plot to {}".format(args.output_stub+".png"))
        plt.savefig(args.output_stub+".png",dpi=200)
        # logging.info("Saving plot to {}".format(args.output_stub+".eps"))
        # plt.savefig(output_stub+".eps")
        return

    # def toggle_errorbar_vis(ebarplot):
    #     for i in flatten(ebarplot):
    #         if i:
    #             i.set_visible(not i.get_visible())

    # def func(label):
    #     toggle_errorbar_vis(plots[label])
    #     if label in tmin_plot.keys():
    #         tmin_plot[label].set_visible(not tmin_plot[label].get_visible())
    #     plt.draw()

    # if not seperate:
    #     rax = plt.axes([0.85, 0.8, 0.1, 0.15])
    #     check = CheckButtons(rax, plots.keys(), [True]*len(plots))
    #     check.on_clicked(func)
    #     leg.draggable()

    plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="plot a set of data files")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="increase output verbosity")
    parser.add_argument("-r", "--real", action="store_true",
                        help="don't include the imgainry part'")
    parser.add_argument("-t", "--title", type=str, required=False,
                        help="plot title", default=None)
    parser.add_argument("-s", "--seperate", action="store_true", required=False,
                        help="plot one column or multi columns")
    parser.add_argument("-y", "--yrange", type=float, required=False, nargs=2,
                        help="set the yrange of the plot", default=None)
    parser.add_argument("-x", "--xrang", type=float, required=False, nargs=2,
                        help="set the xrang of the plot", default=None)
    parser.add_argument("-o", "--output-stub", type=str, required=False,
                        help="stub of name to write output to")
    # parser.add_argument('files', metavar='f', type=argparse.FileType('r'), nargs='+',
    #                     help='files to plot')
    parser.add_argument('files', metavar='f', type=str, nargs='+',
                        help='files to plot')
    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.DEBUG)
        logging.debug("Verbose debuging mode activated")
    else:
        logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.INFO)

    boxplot_files()
