#!/usr/bin/env python2
import matplotlib.pyplot as plt
import matplotlib as mpl
import logging
import argparse
import os
import pandas as pd
import math
from cStringIO import StringIO
import re
import numpy as np
from matplotlib.patches import Ellipse
from matplotlib.patches import Rectangle

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

def get_singles(singlefilename):
    if not os.path.isfile(singlefilename):
        logging.warn("Single file missing, ignoring")
        return []
    else:
        txt = lines_without_comments(singlefilename)
        df = pd.read_csv(txt, sep=",", delimiter=",", skipinitialspace=True, names=["index", "levelnum"])
        return list(df.levelnum.values)


def get_colors(colorfilename):
    if not os.path.isfile(colorfilename):
        logging.warn("color file missing, ignoring")
        return None
    else:
        txt = lines_without_comments(colorfilename)
        df = pd.read_csv(txt, sep=",", delimiter=",", skipinitialspace=True,
                         index_col=0, header=None)
        return (df.transpose()/df.max(axis=1)).max(axis=1)



def read_file(filename):
    with open(filename, 'r') as f:
        first_line = f.readline()
    names = [s.strip(" #") for s in first_line.split(",")[0:-2]]
    txt = lines_without_comments(filename)
    filetype = determine_type(txt)
    if filetype == "paren_complex":
        df = pd.read_csv(txt, delimiter=' ', names=names,
                         converters={1: parse_pair, 2: parse_pair})
    if filetype == "comma":
        df = pd.read_csv(txt, sep=",", delimiter=",", names=names, skipinitialspace=True,
                         delim_whitespace=True, converters={0: removecomma, 1: myconverter})
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

    if all([re.match(".*[0-9]_srcCol[0-9].*", n) for n in names]):
        names = ["level_"+re.match(".*[0-9]+_srcCol([0-9]+).*", n).group(1) for n in names]

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


def format_error_string(value, error):
    digits = -1.0*round(math.log10(error))
    if np.isnan(digits):
        return "****"
    formated_error = int(round(error * (10**(digits + 1))))
    formated_value = "{m:.{d}}".format(d=int(digits) + 1, m=value)
    return "{m}({e})".format(m=formated_value, e=formated_error)


def escale(i):
    if not args.experiment:
        return i
    else:
        return (5.0 * i) / (3.0 * 1.67245)

def tscale(i):
    lattice_omegamass=0.27803
    if not args.experiment:
        return i
    else:
        return (5.0 * i) / (3.0 * lattice_omegamass)


def add_experiment_results(experimental_results, f, ax):

    with open(experimental_results, "r") as expfile:
        for line in expfile:
            if line.startswith("#"):
                continue
            name, mass, uncertainty = [i.strip() for i in  line.split(",")]
            smass = escale(float(mass))
            suncertainty = escale(float(uncertainty))
            loc,tloc = -1.6,-1.8
            if "_3" in name:
                loc,tloc = -1.45, -0.9
            rect = plt.Rectangle((loc, smass-suncertainty), 0.5, 2*suncertainty,
                                 fc='r', fill=True, linewidth=1, color='r')
            f.gca().add_artist(rect)
            ax.annotate("${}$".format(name), xy=(tloc, smass), fontsize=18)
    plt.xlim(-2.0, 1.0)
    plt.ylim(0, 3.0)
    plt.plot([-0.5,-0.5], [0,3], 'k-', lw=2, )
    ax.annotate("Experiment", xy=(-1.5, 0.1), fontsize=18, fontweight='bold')
    ax.annotate("Lattice", xy=(0.2, 0.1), fontsize=18, fontweight='bold')
    ax.set_ylabel("$5m/3m_\omega$", fontweight='bold', fontsize=18)

def boxplot_files():
    markers = ['o', "D", "^", "<", ">", "v", "x", "p", "8"]
    # colors, white sucks
    mpl.rcParams['axes.linewidth'] = 5.0
    colors = [c for c in mpl.colors.colorConverter.colors.keys() if c != 'w' and c != "k"]
    cm = plt.get_cmap("winter_r") # colormap to use
    colors.append("#ffa500")
    plots = {}
    labels = label_names_from_filelist(args.files)
    #labels = [translate(l) for l in labels]
    circles = []
    single_indecies = []
    outline_single = []
    mycolors = []
    data = []
    f, ax = plt.subplots()

    prevtextloc = 0.0
    dfs = {}
    for label, filename in zip(labels, args.files):
        dfs[label] = read_file(filename)

    with open(args.ordering) as orderfile:
        ordering = [i.strip() for i in orderfile.readlines()]

    sdfs = [(i,dfs[i]) for i in ordering]
    sdfs = [i for i in sdfs if i[1].mass.std() < args.prune]
    if args.maxlevels:
        sdfs = sdfs[:args.maxlevels]

    if args.experiment and args.single:
        logging.info("experiment and single, so only plot the singles")

        sdfs = [(l,df) for l,df in sdfs if int(l) in args.single]

    sorted_labels = [i[0] for i in sdfs]

    for index, (label, df) in enumerate(sdfs):

        # for index, label in enumerate(labels):
        color = "b" if args.experiment else colors[index % len(colors)]


        if args.seperate:
            data.append(df.mass.values)
            levelnum=int(label)
            if args.color is not None:
                if args.splitbox:
                    lower = df.mass.quantile(q=0.25)
                    upper = df.mass.quantile(q=0.75)
                    circles.append(Rectangle((index+0.75, lower), width=0.5, height=(upper-lower)*args.color[levelnum+1], color='b', fill=True))
                elif args.outline:
                    if args.color[levelnum+1] == 1.0:
                        single_indecies.append(index)
                    if 1.0 > args.color[levelnum+1] > 0.5:
                        outline_single.append(index)

                else:
                    mycolors.append(args.color[levelnum+1])
            if levelnum in args.single:
                logging.info("adding level{} index {} to single_index".format(levelnum, index))
                single_indecies.append(index)
                circles.append(Ellipse((index+1, df.mass.median()), width=1.1, height=df.mass.std()*5.0, color='r', fill=False))
        else:
            med = tscale(df.mass.median())
            width = tscale(df.mass.std())
            values = tscale(df.mass.values)

            offset = 0.25+((1-(index+1) % 3) * 0.33)#+(index/3)*0.05
            if index%3 == 0 and index%2==0 :
                offset += (index/3)*0.03
            if index%3 == 2 and index%2==0:
                offset += (index/3)*0.03
            if index%3 == 1 and index%2==0:
                offset -= (index/3)*0.03
            prevtextloc = med if med-prevtextloc > 0.01 else prevtextloc+0.01

            textloc = (-1.2 if (index + 1) % 3 > 0 else 1, prevtextloc)
            plots[label] = plt.boxplot(values, widths=0.5, patch_artist=True,
                                       positions=[offset])
            hide = not args.clean
            plots[label]["boxes"][0].set_facecolor(color)
            plots[label]["boxes"][0].set_linewidth(1)
            plots[label]["boxes"][0].set_color(color)
            plots[label]["boxes"][0].set_alpha(max(1.0-width*2.0, 0.1))
            plots[label]["boxes"][0].set_zorder(-1*width)
            plt.setp(plots[label]["whiskers"], color=color, visible=hide)
            plt.setp(plots[label]["fliers"], color=color, visible=hide)
            plt.setp(plots[label]["caps"], color=color, visible=hide)
            plt.setp(plots[label]["medians"], visible=hide)
            if not args.experiment:
                ax.annotate(label+":{}".format(format_error_string(med, width)), xy=(offset-0.1, med),
                            xytext=textloc, arrowprops=dict(arrowstyle="->", fc="0.6"))

    if args.seperate:
        splot = plt.boxplot(data, widths=0.5, patch_artist=True)
        for i,b in enumerate(splot["boxes"]):
            if args.color is not None and not args.splitbox and not args.outline:
                color = cm(mycolors[i]) # colormap
            else:
                color = 'c'
            b.set_linewidth(2)
            b.set_facecolor(color)
            b.set_linewidth(3)
            b.set_color(color)
            # b.set_alpha(mycolors[i])
            if i in single_indecies:
                b.set_facecolor('b')
                b.set_color('b')
            if i in outline_single:
                b.set_color('b')
                b.set_facecolor('c')
        if args.clean:
            plt.setp(splot["whiskers"], visible=False)
            plt.setp(splot["fliers"], visible=False)
            plt.setp(splot["caps"], visible=False)
            plt.setp(splot["medians"], visible=False)
            plt.setp( ax.get_xticklabels(), visible=False)
        else:
            xticknames = plt.setp(ax, xticklabels=sorted_labels)
            plt.setp(xticknames, rotation=45, fontsize=8)
        for c in circles:
            f.gca().add_artist(c)
        ax.set_ylabel("$a_t$Energy", fontweight='bold')
        ax.set_xlabel("Levels", fontweight='bold')
        ax.yaxis.set_tick_params(width=5, length=10)
        ax.xaxis.set_tick_params(width=2, length=6)
    if not args.seperate:
        plt.tick_params(labelbottom="off", bottom='off')
        if args.experiment:
            add_experiment_results(args.experiment, f, ax)
        else:
            plt.xlim(-1.5, 1.5)

    if args.yrange:
        plt.ylim(args.yrange)


    if args.title:
        f.suptitle(args.title.replace("_", " "), fontsize=24)

    if args.threshold:
        plt.plot([-2, 200], [args.threshold, args.threshold], color='r', linestyle='--', linewidth=2)

    if(args.output_stub):
        f.set_size_inches(19.2, 12.0)
        plt.rcParams.update({'font.size': 24})
        f.set_dpi(100)
        #plt.tight_layout(pad=2.0, h_pad=1.0, w_pad=2.0)
        #plt.tight_layout()
        logging.info("Saving plot to {}".format(args.output_stub+".png"))
        plt.savefig(args.output_stub+".png")
        # logging.info("Saving plot to {}".format(args.output_stub+".eps"))
        # plt.savefig(output_stub+".eps")
        return

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
    parser.add_argument("--ordering", type=str, required=True,
                        help="file which contains the ordering")
    parser.add_argument("-o", "--output-stub", type=str, required=False,
                        help="stub of name to write output to")
    parser.add_argument("-y", "--yrange", type=float, required=False, nargs=2,
                        help="set the yrange of the plot", default=None)
    parser.add_argument("-e", "--experiment", type=str, required=False,
                        help="file with experimental results")
    parser.add_argument("-c", "--clean", action="store_true", required=False,
                        help="display without outliers or wiskers")
    parser.add_argument("-3", "--threshold", type=float, required=False,
                        help="Draw a line where 3 particle threshold is")
    parser.add_argument("-n", "--maxlevels", type=int, required=False,
                        help="dont plot more  than this many levels")
    parser.add_argument("-p", "--prune", type=float, required=False, default=1e10,
                        help="remove levels with error above this")
    parser.add_argument("--single", type=str, required=False, default=[],
                        help="file with single hadron info")
    parser.add_argument("--color", type=str, required=False, default=None,
                        help="color code by file")
    parser.add_argument("--splitbox", action="store_true", required=False,
                        help="split the box into color rather than graidant the boxes")
    parser.add_argument("--outline", action="store_true", required=False,
                        help="have color outline")
    # parser.add_argument('files', metavar='f', type=argparse.FileType('r'), nargs='+',
    #                     help='files to plot')
    parser.add_argument('files', metavar='f', type=str, nargs='+',
                        help='files to plot')
    args = parser.parse_args()

    if args.experiment and args.seperate:
        logging.error("Comaparison to experimental results doesn't work in seperate mode'")
        parser.print_help()
        parser.exit()

    if args.single:
        args.single = get_singles(args.single)

    if args.color:
        args.color = get_colors(args.color)

    if args.experiment and args.single:
        logging.info("experiemental and single, so only ploting the singles against experiment")

    if args.verbose:
        logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.DEBUG)
        logging.debug("Verbose debuging mode activated")
    else:
        logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.INFO)

    boxplot_files()
