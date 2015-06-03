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
fitparser.add_argument("-tmax", "--tmax", action="store_true", required=False,
                       help="set end of fit range to max")
fitparser.add_argument("-b", "--bootstraps", type=int, required=False, default=NBOOTSTRAPS,
                    help="Number of straps")
fitparser.add_argument("-p", "--plot", action="store_true", required=False,
                    help="Plot the resulting fit")
fitparser.add_argument("-d", "--debug", action="store_true", required=False,
                       help="debug the fit")
fitparser.add_argument("--histo", action="store_true", required=False,
                       help="make a histogram of the fit paramters")
fitparser.add_argument("--prune", type=float, required=False, default=0.5,
                       help="number of sigma to prune at")
fitparser.add_argument("--symmetric", action="store_true",
                       help="check for symmetry then symmetrize")
fitparser.add_argument("--antisymmetric", action="store_true",
                       help="check for antisymmetry then (anti)symmetrize")
fitparser.add_argument("-Nt", "--period", type=int, required=False,
                    help="Period in time direction (not required for all functions)")
fitparser.add_argument("-r", "--random", type=int, default=4, required=False,
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
fitparser.add_argument("--full", action="store_true",
                       help="if the input correlator is all time slices")
fitparser.add_argument("-m", "--minuit", action="store_true",
                    help="use the minuit fitter")
fitparser.add_argument("--debugguess", action="store_true",
                       help="don't fit, just show the intial guess")
fitparser.add_argument("--debug_identcov", action="store_true",
                       help="set covariance to identity")
fitparser.add_argument("--debug_uncorrelated", action="store_true",
                       help="only perform an uncorrelated fit")
fitparser.add_argument("--jackknife", action="store_true",
                       help="jackknife instead of bootstrap")
fitparser.add_argument("--bin", type=int, required=False,
                       help="bin the correlators first")
fitparser.add_argument("--nofallback", action="store_true",
                       help="don't fall back, for testing")
fitparser.add_argument("-a", "--alltimes", action="store_true", required=False,
                       help="do all possible fits")
fitparser.add_argument("--skip_done", action="store_true", required=False,
                       help="skip the fit if already done")
fitparser.add_argument("-f", "--function", choices=functions.keys(),
                    required=False, default="periodic_exp", help="function to fit to")
