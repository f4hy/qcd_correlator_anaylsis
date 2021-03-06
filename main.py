#!/usr/bin/env python
import plot
import build_corr
import logging
import argparse
import os
import sys
import determine_operators
import fit
import inspect

function_list = inspect.getmembers(sys.modules["fitfunctions"], inspect.isclass)
functions = {name: f for name, f in function_list}

parser = argparse.ArgumentParser(description="compute and plot effective masses")
parser.add_argument("-i", "--input-dir", type=str, required=True,
                    help="directory to read files from")
parser.add_argument("-o", "--output-dir", type=str, required=True,
                    help="directory to write plots to")
parser.add_argument("-ob", "--output-bins", type=str, required=False,
                    help="directory to write binned data to")
parser.add_argument("-r", "--operators", action='append', required=False,
                    help="operator to make \n\n e.g. -r a1pp_0_optype0_op1")
parser.add_argument("-b", "--bins", type=int, default=1, help="number of bins")
parser.add_argument("-m", "--make-from-operators",
                    help="build from operators rather than correlator files", action="store_true")
parser.add_argument("--off-diagonals", action="store_true",
                    help="do the off diagonals")
parser.add_argument("-v", "--verbose", action="store_true",
                    help="increase output verbosity")
parser.add_argument("-f", "--format", type=str, required=False,
                    help="fromat of the correlator files in the directory\n\n"
                    "e.g. {}-{}.A1gp.conn.dat where {} are replaced with operator strings"
                    "Defaults to '{}-{}.dat'",
                    default="{}-{}.dat")
parser.add_argument("-fv", "--format-vev", type=str, required=False,
                    help="fromat of the vev files in the directory\n\n"
                    "e.g. {}-A1gp.conn.vev where {} are replaced with operator strings"
                    "Defaults to '{}.vev'",
                    default=None)
parser.add_argument("-nv", "--no-vev", action="store_true", required=False,
                    help="Specify no vev so should be set to zeros\n")
parser.add_argument("-dt", "--delta-t", nargs='+', required=False, default=[1, 3], type=int,
                    help="which delta-t's to compute for effective masses \n")
parser.add_argument("--function", choices=functions.keys(),
                    required=False, default="single_exp", help="function to fit to (if fitting)")
parser.add_argument("--fit", action="store_true", help="Fit the correaltors, add fit value in the comments")
parser.add_argument("-Nt", "--period", type=int, required=False,
                    help="Period in time direction (not required for all functions)")
parser.add_argument("--periodic", action="store_true", help="use periodic effective mass functions")
parser.add_argument("--prune", action="store_true", help="don't plot stuff where the correlator value is too close to zero (with error)")
parser.add_argument("-c", "--configs", type=int, required=False, help="specify the configs to be used\n")
parser.add_argument("-t", "--times", required=False, type=int, help="specify the times to be used\n")


args = parser.parse_args()

if args.fit:
    raise DeprecationWarning("Fit no long works (plotting fits was reworked)")
    funct = functions[args.function](Nt=args.period)

if not args.operators:
    print "Operators not specified, attempting to automagically determine"
    ops = determine_operators.matching_operators(args.input_dir, args.format)
    print ops
    if not ops:
        print "Error: no operators found"
        parser.print_help()
        parser.exit()
    args.operators = ops

cor_template = args.format
if (not args.make_from_operators) and (not args.format_vev and not args.no_vev):
    print "Error: must specify vev format or no-vev (-nv)"
    parser.print_help()
    print "\nError: must specify vev format or no-vev (-nv)"
    parser.exit()
vev_template = args.format_vev

args.input_dir = os.path.normpath(args.input_dir) + os.sep
args.output_dir = os.path.normpath(args.output_dir) + os.sep

if not os.path.exists(args.input_dir):
    print "input directory doesnt exist"
    parser.print_help()
    parser.exit()


if not args.output_bins:
    args.output_bins = args.output_dir

if args.verbose:
    logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.DEBUG)
    logging.debug("Verbose debuging mode activated")
else:
    logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.INFO)

logging.info("Running with operators" + str([x.strip() for x in args.operators]))


def main():
    if not args.off_diagonals:
        for oper in [unstriped.strip() for unstriped in args.operators]:
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
                correlator = correlator.reduce_to_bins(args.bins)
                correlator.writefullfile(args.output_bins + "binned_%d_%s" % (args.bins, oper))
            if args.fit:
                try:
                    correlator.prune_invalid(delete=True)
                    fitparams = fit.auto_fit(funct, correlator, return_quality=True)
                    plot_corr(correlator, args.output_dir, oper, fitparams)
                except RuntimeError:
                    logging.error("could not fit, skipping fit")
                    plot_corr(correlator, args.output_dir, oper)
            else:
                plot_corr(correlator, args.output_dir, oper)
            logging.info("done with %s %s to %s\n---\n", oper, oper, args.output_dir)
    else:
        for src_oper in args.operators:
            for snk_oper in [unstriped.strip() for unstriped in args.operators]:
                try:
                    if args.make_from_operators:
                        correlator = off_diagonal_ops(args.input_dir, src_oper, snk_oper)
                    else:
                        correlator = off_diagonal_file(args.input_dir, src_oper, snk_oper)

                    if args.bins > 1:
                        binedcor = correlator.reduce_to_bins(args.bins)
                        plot_corr(binedcor, args.output_dir, src_oper + snk_oper)
                        binedcor.writefullfile(args.output_bins + "binned_%d_%s_%s" %
                                               (args.bins, src_oper, snk_oper))
                    else:
                        plot_corr(correlator, args.output_dir, src_oper + snk_oper)
                    logging.info("done with %s %s to %s\n---\n", src_oper, snk_oper, args.output_dir)
                except IOError:
                    logging.error("File not found for {} and {}\nContinuing".format(src_oper, snk_oper))
                    continue


def plot_corr(corr, out_folder, name, fitparams=None):

    if args.prune:
        corr.prune_invalid()

    avgcorr = corr.average_sub_vev()
    corr_errors = corr.jackknifed_errors()

    plot_corr_info = {"%s, \t error" % name: (avgcorr,  corr_errors)}
    plot.plotwitherrorbarsnames("%scorrelator.%s" % (out_folder, name),
                                plot_corr_info, avgcorr.keys(), autoscale=True)

    emass_dts = args.delta_t
    for dt in emass_dts:
        emass = corr.effective_mass(dt)
        emass_errors = corr.effective_mass_errors(dt)
        plot_emass = {"%s emass dt=%d, \t error" % (name, dt): (emass, emass_errors)}
        if fitparams:
            fitcomment = "fit({},{}) m={} e={} qual:{}\n".format(fitparams[0], fitparams[1],
                                                                 fitparams[2][0], fitparams[3][0],
                                                                 fitparams[4])
        else:
            fitcomment = None
        plot.plotwitherrorbarsnames("%semass%d.%s" % (out_folder, dt, name),  plot_emass,
                                    emass.keys(), autoscale=True, addcomment=fitcomment)

        if args.periodic:             # Do it all again with periodic
            cosh_emass = corr.periodic_effective_mass(dt)
            cosh_emass_errors = corr.periodic_effective_mass_errors(dt)
            plot_cosh_emass = {"%s cosh_emass dt=%d, \t error" % (name, dt): (cosh_emass, cosh_emass_errors)}
            plot.plotwitherrorbarsnames("%scosh_emass%d.%s" % (out_folder, dt, name),  plot_cosh_emass,
                                        cosh_emass.keys(), autoscale=True)
            # cosh_const_emass = corr.cosh_const_effective_mass(dt)
            # cosh_const_emass_errors = corr.cosh_const_effective_mass_errors(dt)
            # plot_cosh_const_emass = {"%s cosh_emass dt=%d, \t error" % (name, dt): (cosh_const_emass, cosh_const_emass_errors)}
            # plot.plotwitherrorbarsnames("%scosh_const_emass%d.%s" % (out_folder, dt, name),  plot_cosh_const_emass,
            #                             cosh_const_emass.keys(), autoscale=True)


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
    op_file = data_folder + cor_template.format(op)
    return build_corr.diag_from_opfiles(op_file)


def off_diagonal_ops(data_folder, src_op, snk_op):
    srcop_file = data_folder + cor_template.format(src_op)
    snkop_file = data_folder + cor_template.format(snk_op)
    return build_corr.from_opfiles(srcop_file, snkop_file)


def diagonal_file(data_folder, op):
    corrfile = data_folder + cor_template.format(op, op)
    if(args.no_vev):
        return build_corr.corr_and_vev_from_files(corrfile, cfgs=args.configs, ts=args.times)
    else:
        vev_file = data_folder + vev_template.format(op)
        return build_corr.corr_and_vev_from_files(corrfile, vev_file, vev_file, cfgs=args.configs, ts=args.times)


def off_diagonal_file(data_folder, src_op, snk_op):
    corrfile = data_folder + cor_template.format(src_op, snk_op)
    if(args.no_vev):
        return build_corr.corr_and_vev_from_files(corrfile, cfgs=args.configs, ts=args.times)
    else:
        src_vev_file = data_folder + vev_template.format(src_op)
        snk_vev_file = data_folder + vev_template.format(snk_op)
        return build_corr.corr_and_vev_from_files(corrfile, src_vev_file, snk_vev_file, cfgs=args.configs, ts=args.times)

if __name__ == "__main__":
    main()
