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
import matplotlib.patches as mpatches

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
        logging.debug("paren_complex file type detected")
        return "paren_complex"
    if "," in firstline:
        logging.debug("comma file type detected")
        return "comma"
    logging.debug("space sperated file assumed")
    return "space_seperated"


def get_singles(singlefilename):
    if not os.path.isfile(singlefilename):
        logging.warn("Single file missing, ignoring")
        return []
    else:
        logging.debug("reading {} for levels which are singles".format(singlefilename))
        txt = lines_without_comments(singlefilename)
        df = pd.read_csv(txt, sep=",", delimiter=",", skipinitialspace=True, names=["index", "levelnum"])
        return list(df.levelnum.values)


def get_colors(colorfilename):
    if not os.path.isfile(colorfilename):
        logging.warn("color file missing, ignoring")
        return None
    else:
        logging.debug("reading {} for sh optimize overlaps".format(colorfilename))
        txt = lines_without_comments(colorfilename)
        df = pd.read_csv(txt, sep=",", delimiter=",", skipinitialspace=True,
                         index_col=0, header=None)
        with open(args.ordering) as orderfile:
            ordering = [int(i.strip())+1 for i in orderfile.readlines()]  # ofset by 1 so +1
        todrop = [i for i in df.columns if i not in ordering]
        df = df.drop(todrop, axis=1)
        return (df.transpose()/df.max(axis=1)).max(axis=1)


def read_file(filename):
    with open(filename, 'r') as f:
        first_line = f.readline()
    names = [s.strip(" #") for s in first_line.split(",")[0:-2]]
    if not names:
        names = ["data"]
    txt = lines_without_comments(filename)
    filetype = determine_type(txt)
    if filetype == "paren_complex":
        df = pd.read_csv(txt, delimiter=' ', names=names,
                         converters={1: parse_pair, 2: parse_pair})
    if filetype == "comma":
        df = pd.read_csv(txt, sep=",", delimiter=",", names=names, skipinitialspace=True,
                         delim_whitespace=True, converters={0: removecomma, 1: myconverter, 2:myconverter, 3:myconverter, 4:myconverter, 5:myconverter, 6:myconverter})
    if filetype == "space_seperated":
        df = pd.read_csv(txt, delimiter=' ', names=names)
    return df


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


def format_error_string(value, error):
    digits = -1.0*round(math.log10(error))
    if np.isnan(digits):
        return "****"
    formated_error = int(round(error * (10**(digits + 1))))
    formated_value = "{m:.{d}f}".format(d=int(digits) + 1, m=value)
    return "{m}({e})".format(m=formated_value, e=formated_error)


def scale_params():
    logging.info("setting scale using {}".format(args.scale))
    if args.scale == "none":
        return "$a_t$m", 0.6
    if args.scale == "omega":
        return "$5m/3m_\omega$", 3
    if args.scale == "kaon":
        return "$m/m_K$", 5
    else:
        raise NotImplementedError("Have not implemented {}".format(args.scale))


def experimental_scale(i):
    if args.scale == "none":
        return i
    else:
        if args.scale == "omega":
            return (5.0 * i) / (3.0 * 1.67245)
        if args.scale == "kaon":
            return (i) /  .495
        else:
            raise NotImplementedError("Have not implemented {}".format(args.scale))


def lattice_scale(i):
    if args.scale == "none":
        return i
    if args.scale == "omega":
        lattice_omegamass = 0.27803
        return (5.0 * i) / (3.0 * lattice_omegamass)
    if args.scale == "kaon":
        lattice_kaonmass = 0.0835376129
        return i / lattice_kaonmass
    else:
        raise NotImplementedError("Have not implemented {}".format(args.scale))


def add_experiment_results(experimental_results, f, ax):
    logging.info("Adding experimental results from {}".format(experimental_results))
    with open(experimental_results, "r") as expfile:
        loc, tloc = -2.125, -1.8
        size=0.1
        for line in expfile:
            if line.startswith("#"):
                continue
            name, mass, uncertainty, width, w_uncertainty = [i.strip() for i in line.split(",")]
            smass = experimental_scale(float(mass))
            suncertainty = experimental_scale(float(uncertainty))
            swidth = experimental_scale(float(width))
            sw_uncertainty = experimental_scale(float(w_uncertainty))
            loc+=size+.08
            rect = plt.Rectangle((loc, smass-suncertainty), size, 2*suncertainty,
                                 fc='r', fill=True, linewidth=1, color='r')
            width_rect = plt.Rectangle((loc, smass-swidth), size, 2*(swidth),
                                 fc='r', fill=True, linewidth=1, color='r', alpha=0.3, zorder=-100)
            f.gca().add_artist(rect)
            f.gca().add_artist(width_rect)
            text_yloc= smass-suncertainty-swidth-0.1
            text_yloc= 0.2
            ax.annotate("${}$".format(name), xy=(loc-0.01, text_yloc ), fontsize=21)
    #plt.xlim(-2.0, 1.0)
    label, ymax = scale_params()
    plt.ylim(0, ymax)
    plt.plot([-0.5, -0.5], [0, ymax], 'k-', lw=2, )
    ax.annotate("Experiment", xy=(-1.6, 4.6), fontsize=40, fontweight='bold')
    ax.annotate("Lattice $T_{1u}^+$", xy=(-0.1, 4.6), fontsize=40, fontweight='bold')
    ax.set_ylabel(label, fontweight='bold', fontsize=40)
    return rect, width_rect


def boxplot_files():
    logging.info("makeing boxplots")
    mpl.rcParams['axes.linewidth'] = 5.0
    # colors, white sucks
    colors = [c for c in mpl.colors.colorConverter.colors.keys() if c != 'w' and c != "k"]
    cm = plt.get_cmap("winter_r")  # colormap to use
    colors.append("#ffa500")
    plots = {}
    labels = label_names_from_filelist(args.files)
    circles = []
    single_indecies = []
    outline_single = []
    mycolors = []
    data = []
    f, ax = plt.subplots()

    legend_labels = []

    prevtextloc = 0.0
    dfs = {}
    for label, filename in zip(labels, args.files):
        dfs[label] = read_file(filename)

    if args.ordering:
        with open(args.ordering) as orderfile:
            ordering = [i.strip() for i in orderfile.readlines()]
            logging.debug("Using ordering {}".format(",".join(ordering)))
        sdfs = [(i, dfs[i]) for i in ordering]
    else:
        sdfs = [(i, dfs[i]) for i in labels]

    values = []

    for i in sdfs:
        if len(i[1].columns) == 0:
            continue
        if len(i[1].columns) == 1:
            values.append((i[0], i[1][i[1].columns[0]]))
            continue

        selected_cols = [c for c in i[1].columns if c.startswith(args.col)]
        selected = selected_cols[0]
        if len(selected_cols) > 1:
            logging.warn("More than one {} column, selecint the first one {}".format(args.col, selected_cols))
        #print i[1]["mass"]
        values.append((i[0], i[1][selected]))

    if args.maxlevels:
        sdfs = sdfs[:args.maxlevels]



    sorted_labels = [i[0] for i in sdfs]
    plotindex = -1
    offset = -0.58
    size = 0.125


    for index, (label, value) in enumerate(values):
        # for index, label in enumerate(labels):
        color = "b" if args.experiment else colors[index % len(colors)]
        if args.colorwild:
            color = args.color[index]

        if args.seperate:
            logging.info("Ploting staircase")
            data.append(lattice_scale(value.values))
            try:
                levelnum = int(label)
            except ValueError:
                logging.error("Can not level order, inputs are not levels")
                continue
            if args.color is not None:
                if args.splitbox:
                    lower = lattice_scale(value.quantile(q=0.25))
                    upper = lattice_scale(value.quantile(q=0.75))
                    circles.append(Rectangle((index+0.75, lower), width=0.5,
                                             height=(upper-lower)*args.color[levelnum+1], color='b', fill=True))
                elif args.outline:
                    if args.color[levelnum+1] == 1.0:
                        logging.debug("marking {} as single".format(index))
                        single_indecies.append(index)
                    if 1.0 > args.color[levelnum+1] > args.mixing:
                        logging.debug("marking {} as single-mix".format(index))
                        outline_single.append(index)

                else:
                    mycolors.append(args.color[levelnum+1])
            if levelnum in args.single:
                logging.info("adding level{} index {} to single_index".format(levelnum, index))
                single_indecies.append(index)
                circles.append(Ellipse((index+1, value.median()), width=1.1, height=value.std()*5.0, color='r',
                                       fill=False))
        else:                   # not seperate
            logging.debug("Ploting vertical boxplot")
            if args.experiment and args.color is not None:
                levelnum = int(label)
                if args.color[levelnum+1] < args.mixing:
                    logging.debug("No significant single mixing, skipping")
                    continue
                plotindex += 1
            else:
                plotindex = index

            med = lattice_scale(value.median())
            mean = lattice_scale(value.mean())
            width = lattice_scale(value.std())
            values = lattice_scale(value.values)

            print "level{}: {} {}".format(label, mean, width)
            offset += size+0.1

            prevtextloc = med if med-prevtextloc > 0.01 else prevtextloc+0.01

            textloc = (-1.2 if (plotindex + 1) % 3 > 0 else 1, prevtextloc)
            textloc = (-1.2 if (plotindex + 1) % 3 > 0 else 1, med)
            plots[label] = plt.boxplot(values, widths=size, patch_artist=True,
                                       positions=[offset])
            hide = not args.clean
            b = plots[label]["boxes"][0]
            b.set_facecolor(color)
            b.set_linewidth(1)
            b.set_color(color)
            b.set_alpha(max(1.0-width*2.0, 0.1))
            b.set_zorder(-1*width)
            plt.setp(plots[label]["whiskers"], color=color, visible=hide)
            plt.setp(plots[label]["fliers"], color=color, visible=False)
            plt.setp(plots[label]["caps"], color=color, visible=hide)
            plt.setp(plots[label]["medians"], visible=hide)
            if not args.experiment:
                # ax.annotate(label+":{}".format(format_error_string(med, width)), xy=(offset-0.1, med),
                #             xytext=textloc, arrowprops=dict(arrowstyle="->", fc="0.6"))
                legend_labels.append(mpatches.Patch(color=color, label=label))

            if args.experiment:
                if args.color is not None and 1.0 > args.color[levelnum+1]:
                    logging.debug("marking {} as single-mix".format(index))
                    b.set_linewidth(2)
                    b.set_color('b')
                    b.set_facecolor('c')
                    outlinecolor = b
                else:
                    singlecolor = b
                    b.set_alpha(1)

    if args.seperate:
        splot = plt.boxplot(data, widths=0.5, patch_artist=True)
        for i, b in enumerate(splot["boxes"]):
            if args.color is not None and not args.splitbox and not args.outline:
                color = cm(mycolors[i])  # colormap
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
                singlecolor = b
            if i in outline_single:
                b.set_color('b')
                b.set_facecolor('c')
                outlinecolor = b
            if i not in outline_single and i not in single_indecies:
                normal = b
        if args.clean:
            logging.debug("cleaning up, removing whiskers")
            plt.setp(splot["whiskers"], visible=False)
            plt.setp(splot["fliers"], visible=False)
            plt.setp(splot["caps"], visible=False)
            plt.setp(splot["medians"], visible=False)
            plt.setp(ax.get_xticklabels(), visible=False)
        else:
            xticknames = plt.setp(ax, xticklabels=sorted_labels)
            plt.setp(xticknames, rotation=45, fontsize=8)
        for c in circles:
            f.gca().add_artist(c)
        ax.set_xlabel("Levels", fontweight='bold', fontsize=30)
        ax.xaxis.set_tick_params(width=2, length=6)
    if not args.seperate:
        plt.tick_params(labelbottom="off", bottom='off')
        if args.experiment:
            experimental_box, width_box =add_experiment_results(args.experiment, f, ax)
        else:
            plt.xlim(-1.5, 1.5)
            pass

    if args.ordering:
        ylabel, ymax = scale_params()
        ax.set_ylabel(ylabel, fontweight='bold', fontsize=50)

        if args.yrange:
            plt.ylim(args.yrange)
        else:
            logging.debug("setting yrange to 0,{}".format(ymax))
            # plt.ylim((0, ymax))
        ax.set_yticks(np.arange((ax.get_ylim()[0]), (ax.get_ylim()[1]), 0.2))

    plt.xlim(-1,offset+2)


    if args.outline:
        legend_labels = ["lattice $q \overline{q}$ dominated", "two-hadron dominated", "significant mixing"]
        plt.legend([singlecolor, normal, outlinecolor], legend_labels,
                   fontsize=50, loc=4)
        leg = plt.legend(fancybox=True, shadow=True)
    if args.experiment and args.color is not None:
        # legend_labels = ["single-hadron dominated", "significant mixing", "experimental mass", "experimental width"]
        # plt.legend([singlecolor, outlinecolor, experimental_box, width_box], legend_labels,
        #            fontsize=35, loc=4)
        legend_labels = ["lattice $q \overline{q}$ state", "experimental mass", "experimental width"]
        plt.legend([singlecolor, experimental_box, width_box], legend_labels,
                   fontsize=35, loc=4)
        leg = plt.legend(fancybox=True, shadow=True)
        plt.tick_params(axis='both', which='major', labelsize=50)

    print legend_labels
    leg = plt.legend(handles=legend_labels, loc=0, fontsize=16)

    if args.title:
        f.suptitle(args.title.replace("_", " "), fontsize=50)

    if args.threshold:
        scaled_thresh = experimental_scale(args.threshold)
        logging.info("Drawing line for 3/4 particle threshold at {} {}".format(scaled_thresh, args.threshold))
        plt.plot([-2, 200], [scaled_thresh, scaled_thresh], color='r', linestyle='--', linewidth=2)

    if(args.output_stub):
        f.set_size_inches(20.0, 12.0)
        #plt.rcParams.update({'font.size': 24})
        # plt.tight_layout()
        f.set_dpi(100)
        logging.info("Saving plot to {}".format(args.output_stub+".png"))
        plt.savefig(args.output_stub+".png")
        # ax.set_rasterized(True)
        # logging.info("Saving plot to {}".format(args.output_stub+".eps"))
        # plt.savefig(args.output_stub+".eps")
        # logging.info("Saving plot to {}".format(args.output_stub+".pdf"))
        # plt.savefig(args.output_stub+".pdf")
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
    parser.add_argument("--ordering", type=str, required=False,
                        help="file which contains the ordering")
    parser.add_argument("-o", "--output-stub", type=str, required=False,
                        help="stub of name to write output to")
    parser.add_argument("-y", "--yrange", type=float, required=False, nargs=2,
                        help="set the yrange of the plot", default=None)
    parser.add_argument("-e", "--experiment", type=str, required=False,
                        help="file with experimental results")
    parser.add_argument("--scale",  type=str, required=False, choices=["omega", "kaon", "none"], default="none",
                        help="how to set the sacle")
    parser.add_argument("-c", "--clean", action="store_true", required=False,
                        help="display without outliers or wiskers")
    parser.add_argument("-3", "--threshold", type=float, required=False,
                        help="Draw a line where 3 particle threshold is")
    parser.add_argument("-m", "--mixing", type=float, required=False, default=0.75,
                        help="determine the cutoff to draw the single partile state")
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
    parser.add_argument("--col", type=str, default="mass",
                        help="column of the data file to plot")
    # parser.add_argument('files', metavar='f', type=argparse.FileType('r'), nargs='+',
    #                     help='files to plot')
    parser.add_argument('files', metavar='f', type=str, nargs='+',
                        help='files to plot')
    parser.add_argument("--colorwild", action='append', required=False, default=None,
                        help="color code by wildcard")
    args = parser.parse_args()


    if args.verbose:
        logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.DEBUG)
        logging.debug("Verbose debuging mode activated")
    else:
        logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.INFO)

    if args.experiment and args.seperate:
        logging.error("Comaparison to experimental results doesn't work in seperate mode'")
        parser.print_help()
        parser.exit()

    if args.seperate and args.color:
        if not args.ordering:
            logging.error("seperate color mode requires an ordering")
            parser.print_help()
            parser.exit()

    if args.single:
        args.single = get_singles(args.single)

    if args.color:
        args.color = get_colors(args.color)

    if args.colorwild:
        colors = 'brmckg'
        plotcolors = [""]*len(args.files)
        for c, w in enumerate(args.colorwild):
            for i in range(len(args.files)):
                if w in args.files[i]:
                    plotcolors[i] = colors[c]
        args.color = plotcolors


    if args.experiment and args.color is not None:
        logging.info("experiemental and single, so only ploting the singles against experiment")

    boxplot_files()
