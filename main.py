#!/usr/bin/env python

import plot
import build_corr
import logging
import argparse
import os

parser = argparse.ArgumentParser(description="compute and plot effective masses")
parser.add_argument("-i", "--input-dir", type=str, required=True,
                    help="directory to read files from")
parser.add_argument("-o", "--output-dir", type=str, required=True,
                    help="directory to write plots to")
parser.add_argument("-r", "--operators", action='append', required=True,
                    help="operator to make \n\n e.g. -r a1pp_0_optype0_op1")
parser.add_argument("-b", "--bins", type=int, default=1, help="number of bins")
parser.add_argument("-m", "--make-from-operators",
                    help="build from operators rather than correlator files", action="store_true")
parser.add_argument("--off-diagonals", action="store_true",
                    help="do the off diagonals")
parser.add_argument("-v", "--verbose", action="store_true",
                    help="increase output verbosity")

args = parser.parse_args()

if not os.path.exists(args.input_dir):
    print "input directory doesnt exist"
    parser.print_help()
    parser.exit()

if args.verbose:
    logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.DEBUG)
    logging.debug("Verbose debuging mode activated")
else:
    logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.INFO)

print args.operators


def main():
    if not args.off_diagonals:
        for oper in args.operators:
            if oper == "eigen":
                correlator = eigenvalue_24_balls(args.input_dir)
            elif oper == "eigen32":
                logging.debug("operator eigen32 selected reading 32cubed glueballs")
                correlator = eigenvalue_32_balls(args.input_dir)
            elif args.make_from_operators:
                correlator = diagonal_ops(args.input_dir, oper)
            else:
                correlator = diagonal_file(args.input_dir, oper)

            if args.bins > 1:
                binedcor = correlator.reduce_to_bins(args.bins)
                plot_corr(binedcor, args.output_dir, oper)
                binedcor.writefullfile(args.output_dir + "binned_%d_%s" % (args.bins, oper))
            else:
                plot_corr(correlator, args.output_dir, oper)
            logging.info("done with %s %s to %s\n---\n", oper, oper, args.output_dir)
    else:
        for src_oper in args.operators:
            for snk_oper in args.operators:
                if args.make_from_operators:
                    correlator = off_diagonal_ops(args.input_dir, src_oper, snk_oper)
                else:
                    correlator = off_diagonal_file(args.input_dir, src_oper, snk_oper)

                if args.bins > 1:
                    binedcor = correlator.reduce_to_bins(args.bins)
                    plot_corr(binedcor, args.output_dir, src_oper + snk_oper)
                    binedcor.writefullfile(args.output_dir + "binned_%d_%s_%s" %
                                           (args.bins, src_oper, snk_oper))
                else:
                    plot_corr(correlator, args.output_dir, src_oper + snk_oper)
                logging.info("done with %s %s to %s\n---\n", src_oper, snk_oper, args.output_dir)


def plot_corr(corr, out_folder, name):

    avgcorr = corr.average_sub_vev()
    corr_errors = corr.jackknifed_errors()

    plot_corr_info = {"%s, \t error" % name: (avgcorr,  corr_errors)}
    plot.plotwitherrorbarsnames("%scorrelator.%s" % (out_folder, name),
                                plot_corr_info, avgcorr.keys(), autoscale=True)

    emass_dts = range(1, 4)
    for dt in emass_dts:
        emass = corr.effective_mass(dt)
        emass_errors = corr.effective_mass_errors(dt)
        plot_emass = {"%s emass dt=%d, \t error" % (name, dt): (emass, emass_errors)}

        plot.plotwitherrorbarsnames("%semass%d.%s" % (out_folder, dt, name),  plot_emass,
                                    emass.keys(), autoscale=True)


def eigenvalue_24_balls(data_folder):
    """build a correlator from my 24cubed eigenvalue datafiles
    """
    return build_corr.from_eigenvalue_24cubed_opfiles(data_folder)


def eigenvalue_32_balls(data_folder):
    """build a correlator from my 32cubed eigenvalue datafiles
    """
    return build_corr.from_eigenvalue_32cubed_opfiles(data_folder)


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
    corrfile = "%s%s-%s.A1gp.conn.dat" % (data_folder, op, op)
    vev_file = "%svev.dat" % data_folder
    return build_corr.corr_and_vev_from_files(corrfile, vev_file, vev_file)


def off_diagonal_file(data_folder, src_op, snk_op):
    corrfile = "%s%s-%s.A1gp.conn.dat" % (data_folder, src_op, snk_op)
    src_vev_file = "%svev.dat" % data_folder
    snk_vev_file = "%svev.dat" % data_folder
    # src_vev_file = "%svev_%stest1.dat" % (data_folder, src_op)
    # snk_vev_file = "%svev_%stest1.dat" % (data_folder, snk_op)
    return build_corr.corr_and_vev_from_files(corrfile, src_vev_file, snk_vev_file)

if __name__ == "__main__":
    main()
