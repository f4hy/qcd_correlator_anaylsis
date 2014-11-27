#!/usr/bin/env python
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import logging
import build_corr
import argparse
import fit
from matplotlib.widgets import CheckButtons
from compiler.ast import flatten
import os

from parser_fit import fitparser, functions
import inspect
import sys


NBOOTSTRAPS = 100


def tmin_plot(fn, cor, tmin, tmax, filestub=None, bootstraps=NBOOTSTRAPS):
    emass_dt = 3

    index = fn.parameter_names.index("mass")
    fitted_params = []
    fitted_errors = []
    qualities = []
    Tpoints = range(tmin, tmax-(len(fn.parameter_names)+1))
    if args.write_each_boot:
        orig_write_each_boot = args.write_each_boot
    for t in Tpoints:
        if args.write_each_boot:
            args.write_each_boot = orig_write_each_boot+"_{}".format(t)
        try:
            params, errors, qual = fit.fit(fn, cor, t, tmax,
                                           filestub=filestub, bootstraps=bootstraps, return_quality=True, options=args)
            fitted_params.append(params[index])
            fitted_errors.append(errors[index])
            qualities.append(qual)
        except RuntimeError:
            fitted_params.append(np.nan)
            fitted_errors.append(np.nan)
            qualities.append(0.0)
            continue

    fig = plt.figure()

    emass = cor.effective_mass(emass_dt)
    emass_errors = cor.effective_mass_errors(emass_dt).values()
    emass_plot = plt.errorbar(np.array(emass.keys())+0.2, emass.values(), yerr=emass_errors, fmt='g^', zorder=0)
    cmap = mpl.cm.cool

    tmin_plot = plt.scatter(Tpoints, fitted_params, c=qualities, s=50, cmap=cmap)
    plt.clim(0, 1)
    tmin_error = plt.errorbar(Tpoints, fitted_params, yerr=fitted_errors, fmt=None, zorder=0)

    for i in flatten(emass_plot):
        i.set_visible(False)

    cb = fig.colorbar(tmin_plot)
    cb.set_label("Fit Quality")

    plt.ylim([0, max(emass.values())*1.2])
    plt.xlim([0, tmax + 2])

    def func(label):
        if label == 'tminplot':
            tmin_plot.set_visible(not tmin_plot.get_visible())
            for i in flatten(tmin_error):
                if i:
                    i.set_visible(not i.get_visible())
        elif label == 'emass':
            for i in flatten(emass_plot):
                if i:
                    i.set_visible(not i.get_visible())
        plt.draw()

    if(filestub):
        logging.info("Saving plot to {}".format(filestub+".png"))
        plt.savefig(filestub)
        logging.info("Saving tmin data to {}".format(filestub+".tmin.out"))
        with open(filestub+".tmin.out", "w") as f:
            for t, data, error, q in zip(Tpoints, fitted_params, fitted_errors, qualities):
                f.write("{}, {}, {}, {}\n".format(t, data, error, q))
    else:
        rax = plt.axes([0.85, 0.8, 0.1, 0.15])
        check = CheckButtons(rax, ('tminplot', 'emass'), (True, False))
        check.on_clicked(func)
        plt.show()



if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="compute fits", parents=[fitparser])
    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.DEBUG)
        logging.debug("Verbose debuging mode activated")
    else:
        logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.INFO)

    funct = functions[args.function](Nt=args.period)

    if not args.time_start or not args.time_end:
        parser.error("tmin plot requires specifing a fit range")


    if args.output_stub:
        args.output_stub = os.path.splitext(args.output_stub)[0]

    corrfile = args.inputfile
    vev1 = args.vev
    vev2 = vev1
    if args.vev2:
        vev2 = args.vev2

    try:
        cor = build_corr.corr_and_vev_from_files_pandas(corrfile, vev1, vev2)
    except AttributeError:
        logging.info("Failed to read with pandas, reading normal")
        cor = build_corr.corr_and_vev_from_files(corrfile, vev1, vev2)

    tmin_plot(funct, cor, args.time_start, args.time_end,
              filestub=args.output_stub, bootstraps=args.bootstraps)
