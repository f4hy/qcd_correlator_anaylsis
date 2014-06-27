import argparse
import inspect
import sys
from fitfunctions import *  # noqa

function_list = inspect.getmembers(sys.modules["fitfunctions"], inspect.isclass)
functions = {name: f for name, f in function_list}

NBOOTSTRAPS = 1000

fitparser = argparse.ArgumentParser(add_help=False, description="fit")
fitparser.add_argument("-i", "--inputfile", type=str, required=True,
                    help="Correlator file to read from")
fitparser.add_argument("-o", "--output_stub", type=str, required=False,
                    help="stub of name to write output to")
fitparser.add_argument("-wb", "--write_each_boot", action="store_true", required=False,
                    help="stub of name to write each bootstrap output to")
fitparser.add_argument("-v1", "--vev", type=str, required=False,
                    help="vev file to read from")
fitparser.add_argument("-v2", "--vev2", type=str, required=False,
                    help="vev2 file to read from")
fitparser.add_argument("-ts", "--time-start", type=int, required=False,
                    help="first time slice to start a fit, can be a list of times")
fitparser.add_argument("-te", "--time-end", type=int, required=False,
                    help="last time slice to fit, can be a list of times")
fitparser.add_argument("-max", "--maxrange", action="store_true", required=False,
                    help="fit over the full valid range")
fitparser.add_argument("-b", "--bootstraps", type=int, required=False, default=NBOOTSTRAPS,
                    help="Number of straps")
fitparser.add_argument("-p", "--plot", action="store_true", required=False,
                    help="Plot the resulting fit")
fitparser.add_argument("-Nt", "--period", type=int, required=False,
                    help="Period in time direction (not required for all functions)")
fitparser.add_argument("-r", "--random", type=int, required=False,
                    help="set the random seed")
fitparser.add_argument("-v", "--verbose", action="store_true",
                    help="increase output verbosity")
fitparser.add_argument("--unsafe", action="store_true",
                    help="Code usually exits if something goes wrong "
                    "This option will cause the code to fit anyway.")
fitparser.add_argument("--first_pass", action="store_true",
                    help="Do an uncorrelated chi square first, and use that as the guess")
fitparser.add_argument("--reguess", action="store_true",
                    help="use emass on each bootstrap to set inital guess, otherwise use full ensamble as guess")
fitparser.add_argument("-m", "--minuit", action="store_true",
                    help="use the minuit fitter")
fitparser.add_argument("--debugguess", action="store_true",
                       help="don't fit, just show the intial guess")
fitparser.add_argument("--nofallback", action="store_true",
                       help="don't fall back, for testing")
fitparser.add_argument("-f", "--function", choices=functions.keys(),
                    required=False, default="periodic_exp", help="function to fit to")
