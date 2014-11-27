#!/usr/bin/env python2
import matplotlib.pyplot as plt
import matplotlib as mpl
import logging
import argparse
from matplotlib.widgets import CheckButtons
from compiler.ast import flatten
import re
import os
import pandas as pd
import math
import numpy as np
from scipy import stats
from scipy.optimize import leastsq

from cStringIO import StringIO

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


def iqr(x):
    return stats.scoreatpercentile(x, 75) - stats.scoreatpercentile(x, 25)

def lines_without_comments(filename, comment="#"):
    s = StringIO()
    with open(filename) as f:
        for line in f:
            if not line.startswith(comment):
                s.write(line)
    s.seek(0)
    return s

def removecomma(s):
    return int(s.strip(','))

def myconverter(s):
    print s, np.float(s.strip(','))
    try:
        return np.float(s.strip(','))
    except:
        return np.nan


def read_file(filename):
    with open(filename, 'r') as f:
        first_line = f.readline()
    names = [s.strip(" #") for s in first_line.split(",")[0:-2]]

    txt = lines_without_comments(filename)
    df = pd.read_csv(txt, sep=",", delimiter=",", names=names, skipinitialspace=True,
                         delim_whitespace=True, converters={0: removecomma, 1: myconverter})
    return df


def read_pion_masses(args):
    masses = {}
    for filename in args.files:
        print filename
        df=  read_file(filename)
        mom = int(re.search("mom(\d+)", filename).group(1))
        print mom
        masses[mom] = df.mass
    return masses


def plot_dispersion(args, pion_masses):
    fontsettings = dict(fontweight='bold', fontsize=40)
    markers = ['o', "D", "^", "<", ">", "v", "x", "p", "8"]
    colors = ['b', 'r', 'k', 'm', 'c', 'y']
    plots = {}



    for p,m in pion_masses.iteritems():
        print p, m.mean()

    # def fitfunction(parameter, xdata, ydata):
    #     return

    # original_ensamble_params, success = leastsq(fun, initial_guess, args=(x, y), maxfev=10000)

    XIs = []
    for i in range(len(pion_masses[0])):
        restmass = pion_masses[0][i]
        xvalues = np.array(pion_masses.keys())
        yvalues = np.array([pion_masses[p][i] for p in xvalues])
        print xvalues
        print yvalues
        def fitfun(xi, mom, my):
            # print "FIT FUN"
            # print xi, mom, my
            # print type(xi), type(mom), type(my)
            # print (mom * 2*np.pi / 32 )**2
            return ( restmass**2 + mom*(2*np.pi / 32 )**2 / xi**2  - my**2)
        #fun = lambda xi, mom, my: ( restmass**2 + (mom * 2*np.pi / 32 )**2 / xi**2  - my**2)
        # for x in xvalues:
        #     xi_sq = ((4.0*np.pi**2)  * x) / ((32**2) * (pion_masses[x][i]**2-restmass**2))
        #     print (pion_masses[x][i]**2-restmass**2)
        #     print x, xi_sq, np.sqrt(xi_sq)
        fit_param, success = leastsq(fitfun, 3.44, args=(xvalues, yvalues), maxfev=10000)
        XIs.extend(fit_param)
    print XIs
    print np.mean(XIs), np.std(XIs)
    Xdata = np.linspace(0,max(pion_masses.keys()), 100)
    #exit()
    if args.output_stub:
        with open(args.output_stub+".aniso", "w") as outdata:
            for aniso in XIs:
                outdata.write("{}\n".format(repr(aniso)))


    f, layout = plt.subplots(nrows=1, ncols=1)
    # f.suptitle(r"$ \xi $ = {}".format(np.mean(XIs)), **fontsettings)
    layout.set_xlabel("d$^2$", **fontsettings)
    layout.set_ylabel("$E^2$", **fontsettings)

    dispersionplot = layout
    fitdata = (pion_masses[0].mean()**2 + Xdata*(4.0*np.pi**2) / (32.0**2*np.mean(XIs)**2))
    dispersionplot.plot(Xdata, fitdata, label="fit")
    # for p in layout:
    #     p.tick_params(axis='both', which='major', labelsize=20)

    data = {}
    names = label_names_from_filelist(args.files)
    index = 0
    for p, m in pion_masses.iteritems():
        label = names[index]
        mark = markers[index % len(markers)]
        color = colors[index % len(colors)]
        data[label] = (p,m.median(), iqr(m))
        print data[label]

        print "Wtf"
        # print m
        # print m**2
        msqr = m**2
        print p, m[0] , msqr[0]
        plots[p] = dispersionplot.errorbar(p, msqr.median(), yerr=iqr(msqr), label=label , marker=mark, color=color, ms=10)
        index+=1

    plots["fit"] = dispersionplot.plot()
    if args.yrange:
        plt.ylim(args.yrange)
    if args.xrang:
        plt.xlim(args.xrang)

    plt.tick_params(axis='both', which='major', labelsize=50)

    # if args.title:
    #     f.suptitle(args.title, **fontsettings)



    # for p in layout:
    #     p.legend(fancybox=True, shadow=True, loc=2)

    if(args.output_stub):
        with open(args.output_stub+".dat", "w") as outdata:
            outdata.write("mom, mass, mass_err, \n")
            for k,d in data.iteritems():
                outdata.write(k+",\t"+",".join(map(str,d))+"\n")
        f.set_size_inches(18.5, 10.5)
        plt.rcParams.update({'marker.size': 20})
        #plt.tight_layout(pad=2.0, h_pad=1.0, w_pad=2.0)
        plt.tight_layout()
        if args.eps:
            logging.info("Saving plot to {}".format(args.output_stub+".eps"))
            plt.savefig(args.output_stub+".eps")
        else:
            logging.info("Saving plot to {}".format(args.output_stub+".png"))
            plt.savefig(args.output_stub+".png", dpi=400)
        return


    # def toggle_errorbar_vis(ebarplot):
    #     for i in flatten(ebarplot):
    #         if i:
    #             i.set_visible(not i.get_visible())

    # def func(label):
    #     plots[label].set_visible(not plots[label].get_visible())

    # rax = plt.axes([0.8, 0.3, 0.2, 0.5])
    # check = CheckButtons(rax, plots.keys(), [True]*len(plots))
    # check.on_clicked(func)
    # if not args.nolegend:
    #     leg.draggable()


    plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="plot a set of pion bootstrap fits with anisotropy")
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
    parser.add_argument('files', metavar='f', type=str, nargs='+',
                        help='files to plot')
    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.DEBUG)
        logging.debug("Verbose debuging mode activated")
    else:
        logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.INFO)


    pion_masses = read_pion_masses(args)
    plot_dispersion(args, pion_masses)
