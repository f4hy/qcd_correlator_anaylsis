#!/usr/bin/env python
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
    txt = lines_without_comments(filename)
    filetype = determine_type(txt)
    if filetype == "paren_complex":
        df = pd.read_csv(txt, delimiter=' ', names=["time", "correlator", "error", "quality"],
                         converters={1: parse_pair, 2: parse_pair})
    if filetype == "comma":
        df = pd.read_csv(txt, sep=",", delimiter=",", names=["time", "correlator", "error", "quality"], skipinitialspace=True, delim_whitespace=True, converters={0: removecomma, 1: myconverter})
    if filetype == "space_seperated":
        df = pd.read_csv(txt, delimiter=' ', names=["time", "correlator", "error", "quality"])
    return df


def get_fit(filename, noexcept=False):
    with open(filename) as f:
        for line in f:
            if line.startswith("#fit"):
                logging.info("found fit info: {}".format(line))
                function, tmin, tmax, params, errors = line.split(",")
                fittype = function.split(" ")[0].strip()
                fn = function.split(" ")[1].strip()
                tmin = int(tmin.strip(" (),."))
                tmax = int(tmax.strip(" (),."))
                params = [float(i) for i in params.strip(" []\n").split()  ]
                errors = [float(i) for i in errors.strip(" []\n").split()  ]
                return (fittype, fn, tmin, tmax, params, errors)
    if noexcept:
        return ("single_exp", 0, 1, [0.0,0.0], [0.0,0.0])
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
    funmap = {"two_exp": two_exp, "single_exp": single_exp, "periodic_two_exp": two_exp, "fwd-back-exp": single_exp}
    try:
        fittype, function, tmin, tmax, fitparams, fiterrors = get_fit(filename)
        fun = funmap[function]()
        if fittype == "#fit":
            logging.info("correlator fit info")
            xpoints=np.arange(tmin,tmax,0.3)
            fitpoints = fun.formula(fitparams, xpoints )
            ax.plot(xpoints, fitpoints, ls="dashed", color="r")
            return fun.template.format(*fitparams)
        if fittype == "#fit_emass":
            massindex = fun.parameter_names.index("mass")
            mass= fitparams[massindex]
            masserror = fiterrors[massindex]
            xpoints=np.arange(tmin,tmax,1.0)
            fitpoints = fun.formula(fitparams, xpoints )
            emassfit = []
            dt=3
            for i in range(len(fitpoints))[:-dt]:
                emass = (1.0 / float(dt)) * np.log(fitpoints[i] / fitpoints[i + dt])
                emassfit.append(emass)
            ax.plot(xpoints[:-dt], emassfit, ls="dashed", color="r")
            ax.plot(range(tmin, tmax+1), [mass+masserror]*len(range(tmin, tmax+1)), ls="dashed", color="b")
            ax.plot(range(tmin, tmax+1), [mass-masserror]*len(range(tmin, tmax+1)), ls="dashed", color="b")
            digits = -1.0*round(math.log10(masserror))
            formated_error = int(round(masserror * (10**(digits + 1))))
            formated_mass = "{m:.{d}}".format(d=int(digits) + 1, m=mass)
            return "{m}({e})".format(m=formated_mass, e=formated_error)
    except RuntimeError:
        logging.error("File {} had no fit into".format(filename))


def plot_files(files, output_stub=None, yrange=None, xrang=None, cols=-1, fit=False, real=False, title=None):
    markers = ['o', "D", "^", "<", ">", "v", "x", "p", "8"]
    # colors, white sucks
    colors = [c for c in mpl.colors.colorConverter.colors.keys() if c != 'w' and c != "g"]
    plots = {}
    tmin_plot = {}
    has_colorbar = False
    labels = label_names_from_filelist(files)
    #labels = [translate(l) for l in labels]
    seperate = cols > 0
    layout = None
    if seperate:
        f, layout = plt.subplots(nrows=int(math.ceil(float(len(labels))/cols)), ncols=cols, sharey=True, sharex=True)
    for index, label, filename in zip(range(len(files)), labels, files):
        i = (index)/cols
        j = (index) % cols
        if seperate:
            if len(labels) <= cols:
                axe=layout[j]
            else:
                axe=layout[i][j]

        if fit:
            if seperate:
                fitstring = add_fit_info(filename, ax=axe)
            else:
                fitstring = add_fit_info(filename)
            if fitstring:
                if args.fit_only:
                    print "setting label to {}".format(fitstring)
                    label = fitstring
                else:
                    label += " " + fitstring

        # for index, label in enumerate(labels):
        mark = markers[index % len(markers)]
        color = colors[index % len(colors)]
        df = read_file(filename)
        #print df.head(20)
        #exit()
        time_offset = df.time.values+(index*0.1)
        if seperate:
            time_offset=df.time.values
        logging.debug("%s %s %s", df.time.values, df.correlator.values, df.error.values)
        if any(df["quality"].notnull()):
            logging.info("found 4th column, plotting as quality")
            cmap = mpl.cm.cool
            #print df.correlator.values
            if seperate:
                logging.info("plotting {}  {}, {}".format(label, i, j))
                ax = axe
                ax.set_title(label)
                pass
            else:
                ax = plt
            plots[label] = ax.errorbar(time_offset, df.correlator.values, yerr=df.error.values,
                                       linestyle="none", c=color, marker=mark, label=label,
                                       fmt=None, zorder=0)
            tmin_plot[label] = ax.scatter(time_offset, df.correlator.values, c=df.quality.values,
                                          s=50, cmap=cmap, marker=mark)
            tmin_plot[label].set_clim(0, 1)
            if seperate:
                has_colorbar = True
            if not has_colorbar and not seperate:
                cb = plt.colorbar(tmin_plot[label])  # noqa
                has_colorbar = True
            if yrange:
                plt.ylim(yrange)
            if xrang:
                plt.xlim(xrang)

        else:
            if seperate:
                logging.info("plotting {}  {}, {}".format(label, i, j))
                ax = axe
                ax.set_title(label)
            else:
                ax = plt
            if np.iscomplexobj(df.correlator.values):
                plots[label] = ax.errorbar(time_offset, np.real(df.correlator.values),
                                           yerr=np.real(df.error.values),
                                           linestyle="none", c=color, marker=mark, label=label, ms=10)
                if not real:
                    plots["imag"+label] = ax.errorbar(time_offset, np.imag(df.correlator.values),
                                                      yerr=np.imag(df.error.values),
                                                      markerfacecolor='none',
                                                      linestyle="none", c=color, marker=mark, label=None, ms=10)
            else:
                #print df.correlator.values, df.error.values
                plots[label] = ax.errorbar(time_offset, df.correlator.values, yerr=df.error.values,
                                           linestyle="none", c=color, marker=mark, label=label, ms=12)

            if yrange:
                plt.ylim(yrange)
            else:
                axe.set_ylim(min(df.correlator), max(df.correlator))
            if xrang:
                plt.xlim(xrang)
            else:
                axe.set_xlim(min(df.time), max(df.time))


    if not seperate:
        leg = plt.legend(fancybox=True, shadow=True)
    else:
        plt.tight_layout(pad=0.0, h_pad=0.0, w_pad=0.0)
        if has_colorbar:
            f.subplots_adjust(right=0.95)
            cbar_ax = f.add_axes([0.96, 0.05, 0.01, 0.9])
            f.colorbar(tmin_plot[label], cax=cbar_ax)

    if(output_stub):
        if title:
            f.suptitle(title)
        f.set_size_inches(18.5,10.5)
        plt.rcParams.update({'font.size': 20})
        #plt.tight_layout(pad=2.0, h_pad=1.0, w_pad=2.0)
        plt.tight_layout()
        logging.info("Saving plot to {}".format(output_stub+".png"))
        plt.savefig(output_stub+".png",dpi=100)
        # logging.info("Saving plot to {}".format(output_stub+".eps"))
        # plt.savefig(output_stub+".eps")
        return

    def toggle_errorbar_vis(ebarplot):
        for i in flatten(ebarplot):
            if i:
                i.set_visible(not i.get_visible())

    def func(label):
        toggle_errorbar_vis(plots[label])
        if label in tmin_plot.keys():
            tmin_plot[label].set_visible(not tmin_plot[label].get_visible())
        plt.draw()

    if not seperate:
        rax = plt.axes([0.85, 0.8, 0.1, 0.15])
        check = CheckButtons(rax, plots.keys(), [True]*len(plots))
        check.on_clicked(func)
        leg.draggable()

    plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="plot a set of data files")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="increase output verbosity")
    parser.add_argument("-f", "--include-fit", action="store_true",
                        help="check file for fit into, add it to plots")
    parser.add_argument("-fo", "--fit_only", action="store_true",
                        help="replace_labels with fit info")
    parser.add_argument("-r", "--real", action="store_true",
                        help="don't include the imgainry part'")
    parser.add_argument("-s", "--sort", action="store_true",
                        help="attempt to sort them first")
    parser.add_argument("-c", "--columns", type=int, required=False,
                        help="number of columns to make the plot", default=None)
    parser.add_argument("-n", "--number", type=int, required=False,
                        help="number of correlators to include per plot", default=10000)
    parser.add_argument("-t", "--title", type=str, required=False,
                        help="plot title", default=None)
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


    if args.sort:
        try:
            if args.fit_only:
                fitvalues =[get_fit(i, noexcept=True)[2] for i in args.files]
                s = [x[1] for x in sorted(zip(fitvalues, args.files), key=lambda t: float(t[0]))]
            else:
                s = [x[1] for x in sorted(zip(label_names_from_filelist(args.files), args.files), key=lambda t: int(t[0]))]
            args.files = s
        except Exception as e:
            logging.warn("sorting failed")
            print e
            exit()
        else:
            logging.info("level sorting worked")

    def chunks(l, n):
        """ Yield successive n-sized chunks from l.
        """
        for i in xrange(0, len(l), n):
            yield l[i:i+n]

    for index,chunk in enumerate(chunks(args.files, args.number)):
        ostub=args.output_stub
        if args.output_stub and index > 0:
            ostub="{}_{}".format(args.output_stub, index)
        if args.columns:
            logging.info("Plotting each file as a seperate plot")
            plot_files(chunk, output_stub=ostub,
                       cols=args.columns, yrange=args.yrange, xrang=args.xrang, fit=args.include_fit, real=args.real, title=args.title)
        else:
            plot_files(chunk, output_stub=ostub,
                       yrange=args.yrange, xrang=args.xrang, fit=args.include_fit, real=args.real, title=args.title)
