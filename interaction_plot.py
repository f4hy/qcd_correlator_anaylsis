#!/usr/bin/env python
import matplotlib.pyplot as plt
import logging
import argparse
from matplotlib.widgets import CheckButtons
import glob
import numpy as np
import os

def translate(x):
    namemap = {"A1p": "$A_1^+$", "Ep": "$E^+$", "B1p": "$B_1^+$", "B2p": "$B_2^+$"}
    return [namemap[i] for i in x]


def plot_files(args):
    markers = ['o', "D", "^", "<", ">", "v", "x", "p", "8"]
    # colors, white sucks
    # colors = sorted([c for c in mpl.colors.colorConverter.colors.keys() if c != 'w' and c != "g"])
    colors = ['b', 'r', 'k', 'm', 'c', 'y']
    plots = {}
    tmin_plot = {}
    has_colorbar = False
    fontsettings = dict(fontweight='bold', fontsize=50)
    color = 'k'
    # plt.rcParams.update({'font.size': 20})
    # dirs = glob.glob(args.dataroot+"/psqr*/*_1/")

    channels = {1:("A1p", "Ep"), 2:("A1p", "B1p", "B2p"), 3:("A1p", "Ep"), 4:("A1p", "Ep") }

    for mom in channels.keys():
        plt.figure(mom)
        thisplot = plt.subplot(111)
        biggest_thresh = 0
        for index,irrep in enumerate(channels[mom]):
            mark = markers[index]
            label = irrep
            plotsettings = dict(linestyle="none", marker=mark, label=label, ms=8, elinewidth=3, capsize=20,
                                capthick=2, aa=True, markeredgecolor='none')
            print mom,irrep
            threshold = 0.3
            threshfilename = args.dataroot+"/psqr{}/{}_1/threshold.txt".format(mom,irrep)
            if os.path.isfile(threshfilename):
                with open(threshfilename) as threshfile:
                    threshold = float(threshfile.read())
                logging.info("Adding threshold at {}".format(threshold))
                plt.fill_between([index,index+1], [threshold,threshold], 1.0, color='b', hatch="/", alpha=0.9, linewidth=1)
                biggest_thresh = max(biggest_thresh,threshold)

            noninteractingfilename = args.dataroot+"/psqr{}/{}_1/noninteracting.txt".format(mom,irrep)
            if os.path.isfile(noninteractingfilename):
                with open(noninteractingfilename) as noninteractingfile:
                    for line in noninteractingfile:
                        noninteractinglevel = float(line)
                        logging.info("Adding noninteracting at {}".format(noninteractinglevel))
                        plt.plot([index,index+1], [noninteractinglevel,noninteractinglevel], linestyle="--", color='g', linewidth=4)

            filename = args.dataroot+"/psqr{}/{}_1/diag_fulltops/5_8_cond5000.0/fits/two_exp_tmax.summary".format(mom,irrep)
            with open(filename) as sumfile:
                for line in sumfile:
                    color = colors[index]
                    level, amp, amperr, mass, masserr = map(float,line.split(","))
                    print level
                    if mass < 0 or mass > threshold:
                        continue
                        print mass,masserr
                    if mass+masserr*2 > threshold or (mom == 2 and irrep == "B1p" and level == 1):
                        print "OMG", mass, masserr, threshold
                        color='0.5'
                    plots[label] = plt.errorbar(index+0.5, mass, yerr=masserr, color=color, **plotsettings)
                    plt.plot([index,index], [0,1.0], color='k', lw=4)
        plt.xlim(0,len(channels[mom]))
        plt.xticks(np.array(range(len(channels[mom])))+0.5, **fontsettings)
        thisplot.set_xticklabels(translate(channels[mom]), fontsize=50)
        plt.title("$P^2$={}".format(mom), **fontsettings)
        plt.ylim((0.1,biggest_thresh+0.01))

        plt.ylabel("$a_t E$", fontsize=50)

    # for d in dirs:
    #     if "psqr0" in d:
    #         continue
    #     print d
    #     # print os.path.join(d,"diag_fulltops/5_8_cond5000.0/two_exp_tmax.summary")
    #     filename = os.path.join(d,"diag_fulltops/5_8_cond5000.0/fits/two_exp_tmax.summary")

    # plt.rcParams.update({'font.size': 20})

    if(args.output_stub):
        for mom in channels.keys():
            f = plt.figure(mom)
            f.set_size_inches(18.5, 10.5)
            # plt.rcParams.update({'font.size': 20})
            plt.tight_layout(pad=2.0, h_pad=1.0, w_pad=2.0)
            #plt.tight_layout()
            filename = "{}_{}.png".format(args.output_stub, mom)
            if args.eps:
                filename = "{}_{}.eps".format(args.output_stub, mom)
            logging.info("Saving plot to {}".format(filename))
            plt.savefig(filename)
        return


    plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="plot a set of data files")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="increase output verbosity")
    parser.add_argument("-f", "--include-fit", action="store_true",
                        help="check file for fit into, add it to plots")
    parser.add_argument("-e", "--eps", action="store_true",
                        help="save as eps not png")
    parser.add_argument("-nl", "--nolegend", action="store_true",
                        help="Don't plot the legend")
    parser.add_argument("-fo", "--fit_only", action="store_true",
                        help="replace_labels with fit info")
    parser.add_argument("-ff", "--fitfunction", action="store_true",
                        help="replace_labels with fit function")
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
    parser.add_argument("-tr", "--translate", action="store_true", required=False,
                        help="Attempt to translate the names (of operators)")
    parser.add_argument("-y", "--yrange", type=float, required=False, nargs=2,
                        help="set the yrange of the plot", default=None)
    parser.add_argument("-x", "--xrang", type=float, required=False, nargs=2,
                        help="set the xrang of the plot", default=None)
    parser.add_argument("-o", "--output-stub", type=str, required=False,
                        help="stub of name to write output to")
    # parser.add_argument('files', metavar='f', type=argparse.FileType('r'), nargs='+',
    #                     help='files to plot')
    parser.add_argument("-i", "--dataroot", type=str, required=False,
                        help="root data directory")
    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.DEBUG)
        logging.debug("Verbose debuging mode activated")
    else:
        logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.INFO)

    plot_files(args)
