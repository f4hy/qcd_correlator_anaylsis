#!/usr/bin/env python
import argparse
import logging
import readinput
import os
import re
import irreps
import sys
import particle_operators
from primary_operators import read_coeffs
from prunedops import getprunedops
from momentum_permute import all_permutations

expected_levels_path = "/home/colin/research/notes/hadron_spectrum/expectedlevels/final_results"
#operators_path = "/latticeQCD/raid6/bfahy/operators"
coeffs_path = "/latticeQCD/raid1/laph/qcd_operators/meson_meson_operators/mom_ray_{}"

mom_map = {"#": 2, "+": 1, "0": 0, "-": -1, "=": -2}
cg_map = {"0": "", "1": "CG_1 "}


def rotate(opline):
    print opline
    mom1,mom2 = re.findall("P=\(.*?\)", opline)
    print (mom1, mom2)
    print opline.split()[:8]
    if "CG" in opline:
        base, channel, CG, p1, irrep1, disp1, p2, irrep2, disp2 = opline.split()[:9]
    else:
        base, channel, p1, irrep1, disp1, p2, irrep2, disp2 = opline.split()[:8]
    # print [i for i in mom1.strip("P=()").split(",")]
    # exit()
    m1 = mom1.strip("P=()").split(",")
    m2 = mom2.strip("P=()").split(",")
    psq1 = sum(int(i)**2 for i in m1)
    psq2 = sum(int(i)**2 for i in m2)

    p= "".join([str(int(i)+int(j)) for i,j in zip(m1,m2)])
    print p
    psq_total = sum((int(i)+int(j))**2 for i,j in zip(m1,m2))
    momray = p.replace("-1", "-").replace("-2","=").replace("1","+").replace("2","#")
    print momray
    channel_1 = channel.split("_")[0]+"_1"
    print opline
    coeffs = read_coeffs(momray, channel_1)
    print opline.split()
    expression = ".*{}_.*{}_.*".format(irrep1, irrep2)

    for momentum in all_permutations(psq_total, nobackies=True):
        print momentum
        momray = momentum.replace("-1", "-").replace("-2","=").replace("1","+").replace("2","#")
        print momray
        coeffs = read_coeffs(momray, channel_1)
        print coeffs

        for c in coeffs:
            if re.match(expression, c):
                print c
                S, mom1, i1, mom2, i2, cg = c.split("_")
                coeffpsq1 = sum(mom_map[m]**2 for m in mom1)
                coeffpsq1_2 = sum((mom_map[m]*2)**2 for m in mom1)
                coeffpsq1_3 = sum((mom_map[m]*3)**2 for m in mom1)
                coeffpsq2 = sum(mom_map[m]**2 for m in mom2)
                coeffpsq2_2 = sum((mom_map[m]*2)**2 for m in mom2)
                coeffpsq2_3 = sum((mom_map[m]*3)**2 for m in mom2)
                matches = (coeffpsq1 == int(psq1) and coeffpsq2 == int(psq2))
                scale = 1
                if (coeffpsq1_2 == int(psq1) and coeffpsq2_2 == int(psq2)):
                    matches = True
                    scale = 2
                elif (coeffpsq1_3 == int(psq1) and coeffpsq2_3 == int(psq2)):
                    matches = True
                    scale = 3
                if matches:
                    logging.info("{} matched irreps and moms".format(c))
                    m1 = "({},{},{})".format(*[mom_map[m]*scale for m in mom1])
                    m2 = "({},{},{})".format(*[mom_map[m]*scale for m in mom2])
                    logging.info("Found coeff with psqr{}, psqr{} {}".format(psq1, psq2, c))
                    found = True
                    temp = '{} {} {}[P={} {} {} [P={} {} {}\n'
                    newopline = temp.format(base, channel, cg_map[cg], m1, irrep1, disp1, m2, irrep2, disp2)
                    print newopline
                    print opline
                    with open(args.outputstub+"_{}_{}".format(channel,momray), 'a') as writefile:
                        writefile.write(newopline)



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="rotate operators to other momentum")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="increase output verbosity")
    parser.add_argument("-i", "--inputfile", type=str, required=True,
                           help="file to read operators from")
    parser.add_argument("-o", "--outputstub", type=str, required=True,
                           help="file to write new operators to")

    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(format='# %(levelname)s: %(message)s', level=logging.DEBUG)
        logging.debug("Verbose debuging mode activated")
    else:
        logging.basicConfig(format='# %(levelname)s: %(message)s', level=logging.WARN)

    with open(args.inputfile) as readfile:
        for l in readfile:
            if l.strip().startswith("@"):
                rotate(l)
