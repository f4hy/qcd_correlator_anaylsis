#!/usr/bin/env python
import argparse
import logging
import readinput
import os
import re
import irreps
import sys

expected_levels_path = "/home/colin/research/notes/hadron_spectrum/expectedlevels/final_results"
operators_path = "/latticeQCD/raid6/yjhang/multi_hadron_pruning_operators"
#operators_path = "/latticeQCD/raid6/bfahy/operators"
coeffs_path = "/latticeQCD/raid1/laph/qcd_operators/meson_meson_operators/mom_ray_000"


def component_ops(opline):
    logging.debug("Computing components of {}".format(opline))
    if not "iso" in line:
        logging.debug("line is not a multiop directive")
        return
    types, totalirrep, p1, irrep1, disp1, p2, irrep2, disp2 = line.split()
    _, type1, type2 = types.split("_")
    p1, p2 = p1.strip("[]P="), p2.strip("[]P=")
    disp1, disp2 = disp1.strip("[]\'()"), disp2.strip("[]\"()")
    print '@oplist.push("{} P={} {}_1 {}")'.format(type1, p1, irrep1, disp1)
    print '@oplist.push("{} P={} {}_1 {}")'.format(type2, p2, irrep2, disp2)
    #exit()


def check_lowest(op):
    pass


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="determine the single hadrons we want to know for a multihadron channel")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="increase output verbosity")
    parser.add_argument('oplist', nargs='?', type=argparse.FileType('r'), help="operator list file")

    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(format='# %(levelname)s: %(message)s', level=logging.DEBUG)
        logging.debug("Verbose debuging mode activated")
    else:
        logging.basicConfig(format='# %(levelname)s: %(message)s', level=logging.WARN)


    for line in args.oplist:
        if line.startswith("@"):
            component_ops(line)
