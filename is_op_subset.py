#!/usr/bin/env python
import argparse
import logging
import os
from level_identifier import readops


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="oplist B is a sebset of oplist A")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="increase output verbosity")
    parser.add_argument("-c", "--correlators", type=str, help="location of correlators",
                        default=os.path.join(os.getcwd(), "normed_correlators/"))
    parser.add_argument('A', metavar='oplistA', type=str, help='OplistA')
    parser.add_argument('B', metavar='oplistB', type=str, help='OplistB')
    parser.add_argument("-f", "--format", type=str, required=False,
                        help="fromat of the correlator files in the directory\n\n"
                        "e.g. cor_snk-{}_src-{}.dat where {} are replaced with operator strings"
                        "Defaults to 'ncor_snk-{}_src-{}.dat'",
                        default="ncor_snk-{}_src-{}.dat")

    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(format='# %(levelname)s: %(message)s', level=logging.DEBUG)
        logging.debug("Verbose debuging mode activated")
    else:
        logging.basicConfig(format='# %(levelname)s: %(message)s', level=logging.WARN)

    if set(readops(args.B)) == (set([i for i in readops(args.A) if not i.startswith("iso")])):
        logging.debug("{} is the single particle subset of {}".format(args.B, args.A))
        exit(0)
    else:
        logging.debug("{} is NOT the single subset of {} returning 1".format(args.B, args.A))
        exit(1)
