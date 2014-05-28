#!/usr/bin/env python
import argparse
import logging
import primary_operators
from momentum_permute import all_permutations

expected_levels_path = "/home/colin/research/notes/hadron_spectrum/expectedlevels/final_results"
operators_path = "/latticeQCD/raid6/yjhang/multi_hadron_pruning_operators"
coeffs_path = "/latticeQCD/raid1/laph/qcd_operators/meson_meson_operators/mom_ray_000"

rows = {"A": 1, "B": 1, "E": 2, "T": 3}

def component_ops(opline):
    logging.debug("Computing components of {}".format(opline))
    if "iso" not in line:
        logging.debug("line is not a multiop directive")
        return
    types, totalirrep, p1, irrep1, disp1, p2, irrep2, disp2 = line.split()
    _, type1, type2 = types.split("_")
    p1, p2 = p1.strip("[]P="), p2.strip("[]P=")
    disp1, disp2 = disp1.strip("[]\'()"), disp2.strip("[]\"()")
    if type1 == "kbar":
        type1 = "kaon"
    if type2 == "kbar":
        type2 = "kaon"
    comment1 = "# Lowest state in channel not expected to be a single hadron"
    comment2 = comment1
    def print_op(ptype, p, irrep, disp):
        try:
            lowest_single_index = check_lowest(ptype, p, irrep)
            if lowest_single_index:
                comment = " # There are {} states below the single hadron in this channel".format(lowest_single_index)
            else:
                comment = ""
            if not duplicate(ptype, p, irrep):
                for row in range(1,rows[irrep[0]]+1):
                    print '@oplist.push("{} P={} {}_{} {}"){}'.format(ptype, p, irrep, row, disp, comment)
        except NotImplementedError as e:
            logging.error("{}".format(e))

    psq1 = compute_psq(p1)
    psq2 = compute_psq(p2)
    for m1 in all_permutations(psq1, outputformat="({},{},{})"):
        print_op(type1, m1, irrep1, disp1)
    for m2 in all_permutations(psq2, outputformat="({},{},{})"):
        print_op(type2, m2, irrep2, disp2)


def infer_strangeness_and_isospin(optype):
    namemap = {"kaon": (1, "1h"), "pion": (0, 1), "eta": (0, 0), "phi": (0, 0)}
    return namemap[optype]


def mommap(mom):
    zeroes = mom.count("0")
    ones = mom.count("1")
    twos = mom.count("2")
    if zeroes == 3:
        return "000"
    if ones == 3:
        return "111"
    if ones == 2 and zeroes == 1:
        return "011"
    if ones == 1 and zeroes == 2:
        return "001"
    if twos == 1 and zeroes == 2:
        return "002"
    raise NotImplementedError("unsupported momentum {}".format(mom))


def check_lowest(optype, mom, irrep):
    logging.debug("checking if lowest in channel for {} {} {}".format(optype, mom, irrep))

    strangeness, isospin = infer_strangeness_and_isospin(optype)
    momstring = mommap(mom)
    elevels = primary_operators.read_expected_levels(strangeness, isospin, irrep, mom=momstring)
    index_first_single = next(i for (i,e) in enumerate(elevels) if "PSQ" not in e)
    return index_first_single

duplist = []


def compute_psq(mom):
    nums = mom.strip("()").split(",")
    return sum((int(i)**2 for i in nums))


def duplicate(optype, mom, irrep):
    o = (optype, mom, irrep)
    if o in duplist:
        logging.warning("{} already used!".format(o))
        return True
    else:
        duplist.append((optype, mom, irrep))
        return False


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="determine the single hadrons we want"
                                     "to know for a multihadron channel")
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
        if line.startswith("@") and "CG_1" not in line:
            component_ops(line)
