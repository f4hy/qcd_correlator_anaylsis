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

mom_map = {"#": 2, "+": 1, "0": 0, "-": -1, "=": -2}


def not_implemented(description, default=""):
    logging.error("No primary operator set for particle {}!".format(description))
    logging.error("Unable to handle particles which have"
                  " no primary operator set when run in script mode")
    return None

readinput.askoperator = not_implemented


def single_hadrons():
    logging.info("Include all the single hadrons")
    args.outfile.write("# start with all the single hadrons\n") #beyonce
    if args.isospin == "1":
        singlepath = os.path.join(args.opsdir, "pion/SH")
        with open(singlepath, 'r') as singlefile:
            for line in singlefile:
                args.outfile.write(line)
    if args.isospin == "1h":
        singlepath = os.path.join(args.opsdir, "kaon/SH")
        with open(singlepath, 'r') as singlefile:
            for line in singlefile:
                args.outfile.write(line)
    if args.isospin == "0":     # if I=0 we need eta ops and phi ops
        singlepath = os.path.join(args.opsdir, "eta/SH")
        with open(singlepath, 'r') as singlefile:
            for line in singlefile:
                args.outfile.write(line)
        singlepath = os.path.join(args.opsdir, "phi/SH")
        with open(singlepath, 'r') as singlefile:
            for line in singlefile:
                args.outfile.write(line)

def custom_psqlevel(level, psqr, p1, p2, p1flavor, p2flavor, channel, outfile):
    cg_map = {"0": "", "1": "CG_1 "}
    isomap = {"0": "isosinglet", "1": "isotriplet", "1h": "isodoublet"}
    isoterm = isomap[args.isospin]
    coeffsdir = os.path.join(coeffs_path, channel)
    logging.info("Looking for coeffs in {}".format(coeffsdir))
    coeffs = os.listdir(coeffsdir)
    expression = ".*{}_.*{}_.*|.*{}_.*{}_.*".format(p1[2], p2[2], p2[2], p1[2])
    found = False
    for c in coeffs:
        if re.match(expression, c):
            mom1 = c.split("_")[1]
            mom2 = c.split("_")[3]
            cg = c.split("_")[-1]
            if sum(mom_map[m]**2 for m in mom1) == int(psqr):
                m1 = "({},{},{})".format(*[mom_map[m] for m in mom1])
                m2 = "({},{},{})".format(*[mom_map[m] for m in mom2])
                logging.info("Found coeff with psqr{} {}".format(psqr, c))
                found = True
                temp = '@oplist.push("{}_{}_{} {} {}[P={} {} SS_0] [P={} {} SS_0]")\n'
                opline = temp.format(isoterm, p1flavor, p2flavor, channel, cg_map[cg], m1, p1[2], m2, p2[2])
                args.outfile.write(opline)

    if not found:
        logging.critical("Could not find psqr{} coeffs for this level".format(psqr))
        args.outfile.write("# Could not find psqr{} coeffs for this level".format(psqr))


def flavor_type(particle):
    if particle.I is 1:
        return "pion"
    if particle.I is 0:
        return "eta"
    if particle.I == "1/2":
        return "kaon"


def folder_from_flavor(particle1, particle2):
    pass


def swap(particle1, particle2):
    order = ["kaon", "eta", "phi", "pion"]
    p1 = flavor_type(particle1)
    p2 = flavor_type(particle2)
    if order.index(p1) > order.index(p2):
        return True
    if p1 == p2 == "kaon":
        return True
    return False


def get_unspecified_parameters(args):

    if not args.isospin:
        print("Select isospin")
        args.isospin = readinput.selectchoices(["0", "1", "1h"], default="1")

    if not args.strangeness:
        print("Select strangeness")
        args.strangeness = readinput.selectchoices(["0", "1", "2"], default="0")

    channeldir = os.path.join(operators_path, "BI{I}S{S}".format(I=args.isospin,
                                                                 S=args.strangeness))
    logging.debug("channel dir is {}".format(channeldir))
    channel_list = os.listdir(channeldir)
    if args.channel:
        if args.channel not in channel_list:
            logging.critical("format of input channel is not correct!"
                             " use e.g. {}".format(channel_list[0]))
            parser.exit()
    else:
        print("Select Channel")
        args.channel = readinput.selectchoices(sorted(channel_list))

    args.opsdir = os.path.join(channeldir, args.channel)


def read_expected_levels(args):
    if args.thirtytwo:
        basedir = os.path.join(expected_levels_path, "32^3_phys/mom_000")
    else:
        basedir = os.path.join(expected_levels_path, "24^3_phys/mom_000")

    if args.isospin == "1h":
        filename = "bosonic_2I=1_S={}_levels.txt".format(args.strangeness)
    else:
        filename = "bosonic_I={}_S={}_levels.txt".format(args.isospin, args.strangeness)

    filepath = os.path.join(basedir, filename)
    logging.info("opening {}".format(filepath))
    expected_level_file = open(filepath, "r")
    chan = args.channel.split("_")[0]
    start = False
    expectedleveltxt = ""
    for line in expected_level_file:
        if chan in line:
            start = True
        if start:
            expectedleveltxt += line
            if not line.strip():
                break

    expected_levels = [level.strip() for level in re.findall("\) (.*)", expectedleveltxt)]
    return expected_levels


def get_ops(args, expected_levels):
    level_num = 1
    already_added = []
    for level in expected_levels:
        if args.review:
            raw_input("Look OK?")
        args.outfile.write("\n# level {} {} \n".format(level_num, level))
        level_num += 1
        try:
            opset = irreps.translate_name_to_irrep(level)
        except NotImplementedError, e:
            logging.warn("This level {} is not supported {}, skipping".format(level, e))
            args.outfile.write("# level is not supported: {}, skipping \n".format(e))
            continue
        for op in opset:
            p1, p2 = op
            if p1[3] is None or p2[3] is None:
                args.outfile.write("# {} contains particle which primary operator"
                                   " was never defined! \n".format(level))
                continue
            if swap(p1[0], p2[0]):
                logging.info("swapping")
                p1, p2 = p2, p1
            flavor1 = flavor_type(p1[0])
            flavor2 = flavor_type(p2[0])
            if flavor2 == "kaon":
                flavor2 = "kbar"
            if (flavor1,p1[2]) == ("eta","A1gp") or (flavor2,p2[2]) == ("eta","A1gp"):
                print "OZI SUPRESSED"
                args.outfile.write("# {} SHOULD BE OZI SUPRESSED!! \n".format(level))
                continue
            flavor_folder = "_".join((flavor1, flavor2))
            opdir = os.path.join(args.opsdir, flavor_folder)
            if "PSQ5" in level:
                logging.warn("PSQ5 level, running custom write")
                custom_psqlevel(level, 5, p1, p2, flavor1, flavor2, args.channel, args.outfile)
                continue
            if "PSQ6" in level:
                logging.warn("PSQ6 level, running custom write")
                custom_psqlevel(level, 6, p1, p2, flavor1, flavor2, args.channel, args.outfile)
                continue
            filename = "S={}_{}_{}_{}_{}_0".format(args.strangeness, p1[1], p1[2], p2[1], p2[2])
            filepath = os.path.join(opdir, filename)
            logging.info("opening {}".format(filepath))
            try:
                opfile = open(filepath, 'r')
                opexpression = ".*{} {}.*{} {}.*".format(p1[2], p1[3], p2[2], p2[3])
                momexpression = ".*"
                if "PSQ1" in level or "PSQ2" in level:
                    momexpression = ".*(,1|1,).*"
                if "PSQ4" in level:
                    momexpression = ".*(,2|2,).*"
                if "PSQ9" in level:
                    momexpression = ".*(,3|3,).*"
                found_something = False
                for line in opfile:
                    if re.match(opexpression, line) and re.match(momexpression, line):
                        found_something = True
                        if line in already_added:
                            logging.warn("This operator already added!")
                            args.outfile.write("# This operator already added! \n")
                            args.outfile.write("# {}".format(line))
                        else:
                            args.outfile.write(line)
                            already_added.append(line)
                        if "eta" in line:
                            args.outfile.write("# put in ss operators for every uu operator\n")
                            philine = line.replace("eta", "phi")
                            if philine in already_added:
                                logging.warn("This operator already added!")
                                args.outfile.write("# This operator already added! \n")
                                args.outfile.write("# {}".format(philine))
                            else:
                                args.outfile.write(philine)
                                already_added.append(philine)
                if not found_something:
                    logging.critical("Found no operators for this level!!")
                    logging.critical("Quitting prematurely")
                    exit(-1)
            except IOError:
                logging.info("{} didn't exist".format(filename))
                args.outfile.write("# {} didn't exist\n".format(filename))
    return already_added


def secondary(ops, count):
    args.outfile.write("\n# SECONDARY OPERATORS----------------\n\n")
    newops = []
    for op in ops[:count]:
        words = op.split()
        base, chan, mom1, irrep1, op1, mom2, irrep2, op2 = words
        op1, op2 = op1.strip(']")'), op2.strip(']")')
        _, type1, type2 = base.split("_")
        d = os.path.join(args.opsdir, "_".join((type1, type2)))
        psqr_1 = sum([int(i)**2 for i in mom1[4:-1].split(",")])
        psqr_2 = sum([int(i)**2 for i in mom2[4:-1].split(",")])
        m1 = irreps.momentums[psqr_1]
        m2 = irreps.momentums[psqr_2]
        filename = "S={}_{}_{}_{}_{}_0".format(args.strangeness, m1, irrep1, m2, irrep2)
        #filename = "S={}_{}_{}_{}_{}_0".format(args.strangeness, p1[1], p1[2], p2[1], p2[2])
        filepath = os.path.join(d, filename)
        irrepexpression = ".*{}.*{}.*".format(irrep1, irrep2)
        op1_choices = set()
        op2_choices = set()
        with open(filepath, "r") as f:
            for line in f:
                if re.match(irrepexpression, line):
                    split = line.split()
                    newop1, newop2 = split[4], split[7]
                    op1_choices.add(newop1.strip("]"))
                    op2_choices.add(newop2.strip(']")'))
        # op1_choices.remove(op1)
        # op2_choices.remove(op2)
        print "op1 choices {}".format(" ".join(op1_choices))
        print "op2 choices {}".format(" ".join(op2_choices))
        while True:
            print "Select secondary operator for this level, the primary is:\n {}".format(op)
            print "newop1 {} {} {}? primary is {}".format(type1, mom1, irrep1, op1)
            secondaryop1 = readinput.selectchoices(op1_choices, default=op1)
            print "newop2 {} {} {}? primary is {}".format(type2, mom2, irrep2, op2)
            secondaryop2 = readinput.selectchoices(op2_choices, default=op2)
            secondaryopline = " ".join((base, chan, mom1, irrep1, secondaryop1+"]",
                                        mom2, irrep2, secondaryop2+']")'))+"\n"
            if secondaryopline in ops or secondaryopline in newops:
                print("That choice of ops, already exists in the primary set!!"
                      " Pick a different combination")
            elif secondaryopline in newops:
                print("That choice of ops, already exists in the secondary set!!"
                      " Pick a different combination")
            else:
                logging.info("Secondary operator accepted!")
                break
        newops.append(secondaryopline)
        args.outfile.write(secondaryopline)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="compute and plot effective masses")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="increase output verbosity")
    parser.add_argument("-r", "--review", action="store_true",
                        help="review each operator")
    parser.add_argument("-I", "--isospin", choices=["0", "1", "1h"],
                        help="select isospin")
    parser.add_argument("-S", "--strangeness", choices=["0", "1", "2"],
                        help="select strangeness")
    parser.add_argument("-C", "--channel", type=str,
                        help="select channel e.g. A1up_1")
    parser.add_argument('outfile', nargs='?', type=argparse.FileType('w'),
                        default=sys.stdout)
    parser.add_argument("-32", "--thirtytwo", action="store_true",
                        default=False, help="use 32^3 expected levels")
    parser.add_argument("--secondary", type=int, default=0,
                        help="Prompt user for seconday operators,"
                        " the number is how many operators to make secondaries for")

    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(format='# %(levelname)s: %(message)s', level=logging.DEBUG)
        logging.debug("Verbose debuging mode activated")
    else:
        logging.basicConfig(format='# %(levelname)s: %(message)s', level=logging.WARN)

    get_unspecified_parameters(args)

    single_hadrons()

    expected_levels = read_expected_levels(args)
    ops = get_ops(args, expected_levels)
    if args.secondary:
        secondary(ops, args.secondary)
