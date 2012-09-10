#!/usr/bin/env python

import plot
import build_corr
import correlator
import logging

logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.DEBUG)

BIN_SIZE = 100

ops = ["a1pp_0_optype0_op1", "a1pp_0_optype10_op1"]

#folders
f_data = "/example/"
f_output = "/with/slashes/"


def plot_corr(corr, out_folder, name):


    avgcorr = corr.average_sub_vev()
    corr_errors = corr.jackknifed_errors()

    plot_corr_info = {"%s, \t error" % name: (avgcorr,  corr_errors)}
    plot.plotwitherrorbarsnames("%scorrelator%s" % (out_folder, name),
                                plot_corr_info, avgcorr.keys(), autoscale=True)

    emass_dts = range(1, 4)
    for dt in emass_dts:
        emass = corr.effective_mass(dt)
        emass_errors = corr.effective_mass_errors(dt)
        plot_emass = {"%s emass dt=%d, \t error" % (name, dt): (emass, emass_errors)}

        plot.plotwitherrorbarsnames("%semass%d.%s" % (out_folder, dt, name),  plot_emass,
                                    emass.keys(), autoscale=True)


def diagonal_ops(data_folder, op):
    """return correlator
    """
    op_file = "%sopvals_%stest1.dat" % (data_folder, op)
    return build_corr.diag_from_opfiles(op_file)


def off_diagonal_ops(data_folder, src_op, snk_op):
    srcop_file = "%sopvals_%stest1.dat" % (data_folder, src_op)
    snkop_file = "%sopvals_%stest1.dat" % (data_folder, snk_op)
    return build_corr.from_opfiles(srcop_file, snkop_file)


def diagonal_file(data_folder, op):
    corrfile = "%scor_src_%s-snk_%stest1.dat" % (data_folder, op, op)
    vev_file = "%svev_%stest1.dat" % (data_folder, op)
    return build_corr.corr_and_vev_from_files(corrfile, vev_file, vev_file)

for oper in ops:
    correlator = diagonal_file(f_data, oper)
    #correlator = diagonal_ops(f_data, oper)
    binedcor = correlator.reduce_to_bins(BIN_SIZE)
    plot_corr(binedcor, f_output, oper)
    logging.info("done with %s %s to %s\n\n", oper, oper, f_output)
