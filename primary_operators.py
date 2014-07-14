#!/usr/bin/env python
import argparse
import logging
import readinput
import os
import re
import irreps
import sys
import particle_operators
from prunedops import getprunedops
from momentum_permute import all_permutations

expected_levels_path = "/home/colin/research/notes/hadron_spectrum/expectedlevels/final_results"
#operators_path = "/latticeQCD/raid6/bfahy/operators"
coeffs_path = "/latticeQCD/raid1/laph/qcd_operators/meson_{}_operators/mom_ray_{}"

mom_map = {"#": 2, "+": 1, "0": 0, "-": -1, "=": -2}


def not_implemented(description, default=""):
    logging.error("No primary operator set for particle {}!".format(description))
    logging.error("Unable to handle particles which have"
                  " no primary operator set when run in script mode")
    return None

readinput.askoperator = not_implemented

def single_hadron(ops):
    p = "({},{},{})".format(*[mom_map[i] for i in args.momray])
    isomap = {"1": ["pion"], "1h": ["kaon"], "0": ["eta", "phi"]}
    disp = ops[args.channel.split("_")[0]]
    for name in isomap[args.isospin]:
        singleopline = '@oplist.push("{} P={} {} {}")\n'.format(name, p, args.channel, disp)
        args.outfile.write(singleopline)

def all_single_hadrons(mom):
    logging.info("Include all the single hadrons")
    args.outfile.write("# start with all the single hadrons\n") #beyonce
    isomap = {"1": ["pion"], "1h": ["kaon"], "0": ["eta", "phi"]}

    basechan = args.channel.split("_")[0]

    psqr = sum(int(i)**2 for i in args.momentum)
    if psqr == 0:
        sh = getprunedops(args.isospin, "sh")
        atrest = getprunedops(args.isospin, "atrest")
        ops = sh[basechan] + atrest[basechan]
    if psqr == 1 or psqr == 4:
        ops = getprunedops(args.isospin, "OA")[basechan]
    if psqr == 2:
        ops = getprunedops(args.isospin, "PD")[basechan]
    if psqr == 3:
        ops = getprunedops(args.isospin, "CD")[basechan]
    isomap = {"1": ["pion"], "1h": ["kaon"], "0": ["eta", "phi"]}
    namemap = {"1": "pi", "1h": "K", "0": "eta"}
    psqr = sum(int(i)**2 for i in args.momentum)
    operators = particle_operators.particleDatabase()
    op = operators.read_op(namemap[args.isospin], basechan, psqr)
    for disp in ops:

        for name in isomap[args.isospin]:
            p = "({},{},{})".format(*[mom_map[i] for i in mom])
            singleopline = '@oplist.push("{} P={} {} {}")\n'.format(name, p, args.channel, disp)
            args.outfile.write(singleopline)


def custom_psqlevel(level, psqr, p1, p2, p1flavor, p2flavor, channel, outfile):
    cg_map = {"0": "", "1": "CG_1 "}
    isomap = {"0": "isosinglet", "1": "isotriplet", "1h": "isodoublet", "2": "isoquintet"}
    isoterm = isomap[args.isospin]
    coeffsdir = os.path.join(coeffs_path.format(args.hadrontype, args.momray), channel)
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
        args.outfile.write("# Could not find psqr{} coeffs for this level {}]\n".format(psqr,level))


def read_coeffs(momray, channel):
    coeffsdir = os.path.join(coeffs_path.format(args.hadrontype, momray), channel)
    logging.info("Looking for coeffs in {}".format(coeffsdir))
    coeffs = os.listdir(coeffsdir)
    logging.info("Found {} coeffs".format(len(coeffs)))
    return coeffs


def custom_opline(level, psqr1, psqr2, p1, p2, p1flavor, p2flavor, channel, outfile):
    cg_map = {"0": "", "1": "CG_1 "}
    isomap = {"0": "isosinglet", "1": "isotriplet", "1h": "isodoublet", "2": "isoquintet"}
    isoterm = isomap[args.isospin]

    coeffs = read_coeffs(args.momray, channel)
    if "#" in args.momray:
        down_scaled_coeffs = read_coeffs(args.momray.replace("#","+"), channel)
        coeffs += down_scaled_coeffs
    if p1flavor[0] == p2flavor[0]: # just check fist char for kaon == kbar
        expression = ".*{}_.*{}_.*|.*{}_.*{}_.*".format(p1[2], p2[2], p2[2], p1[2])
    else:
        expression = ".*{}_.*{}_.*".format(p1[2], p2[2], p2[2], p1[2])
    logging.info("searching for {}".format(expression))
    found = False
    oplines = []
    for c in coeffs:
        if re.match(expression, c):
            S, mom1, i1, mom2, i2, cg = c.split("_")
            coeffpsq1 = sum(mom_map[m]**2 for m in mom1)
            coeffpsq1_2 = sum((mom_map[m]*2)**2 for m in mom1)
            coeffpsq1_3 = sum((mom_map[m]*3)**2 for m in mom1)
            coeffpsq2 = sum(mom_map[m]**2 for m in mom2)
            coeffpsq2_2 = sum((mom_map[m]*2)**2 for m in mom2)
            coeffpsq2_3 = sum((mom_map[m]*3)**2 for m in mom2)
            matches = (coeffpsq1 == int(psqr1) and coeffpsq2 == int(psqr2))
            check_swapped = p1flavor[0] == p2flavor[0]
            check_swapped = False
            swapped = False
            scale = 1
            if (coeffpsq1_2 == int(psqr1) and coeffpsq2_2 == int(psqr2)):
                matches = True
                scale = 2
            elif (coeffpsq1_3 == int(psqr1) and coeffpsq2_3 == int(psqr2)):
                matches = True
                scale = 3
            elif check_swapped and (coeffpsq1 == int(psqr2) and coeffpsq2 == int(psqr1)):
                matches = True
                swapped = True
            elif check_swapped and (coeffpsq1_2 == int(psqr2) and coeffpsq2_2 == int(psqr1)):
                matches = True
                swapped = True
                scale = 2
            elif check_swapped and (coeffpsq1_3 == int(psqr2) and coeffpsq2_2 == int(psqr1)):
                matches = True
                swapped = True
                scale = 3

            if matches:
                logging.info("{} matched irreps and moms".format(c))
                if swapped:
                    logging.warn("SWITCHING BECAUSE SWAPMATCH !!!!!!!!!!!!!!!!!!!!!!!!!1")
                    p1, p2 = p2, p1
                m1 = "({},{},{})".format(*[mom_map[m]*scale for m in mom1])
                m2 = "({},{},{})".format(*[mom_map[m]*scale for m in mom2])
                logging.info("Found coeff with psqr{}, psqr{} {}".format(psqr1, psqr2, c))
                found = True
                temp = '@oplist.push("{}_{}_{} {} {}[P={} {} {}] [P={} {} {}]")\n'
                opline = temp.format(isoterm, p1flavor, p2flavor, channel, cg_map[cg], m1, p1[2], p1[3], m2, p2[2], p2[3])
                oplines.append(opline)
            else:
                logging.info("{} matched irreps but not mom".format(c))

    if not found:
        logging.critical("Could not find psqr{}, psqr{} coeffs for this level\n".format(psqr1, psqr2))
        args.outfile.write("# Could not find psqr{},{} ,psqr{},{} coeffs for this level {}]\n".format(psqr1, p1[2], psqr2, p2[2],level))
        return None
    else:
        return oplines

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

    if not args.baryon:
        print("Select hadron type")
        args.hadrontype = readinput.selectchoices(["meson", "baryon"], default="meson")
    else:
        args.hadrontype = "baryon"

    if not args.isospin:
        print("Select isospin")
        args.isospin = readinput.selectchoices(["0", "1", "1h", "2"], default="1")

    if not args.strangeness:
        print("Select strangeness")
        args.strangeness = readinput.selectchoices(["0", "1", "2"], default="0")


    channel_list = os.listdir(coeffs_path.format(args.hadrontype, args.momray))
    if args.channel:
        if args.channel not in channel_list:
            logging.critical("format of input channel is not correct!"
                             " use e.g. {}".format(channel_list[0]))
            parser.exit()
    else:
        print("Select Channel")
        args.channel = readinput.selectchoices(sorted(channel_list))


def read_expected_levels(strangeness, isospin, channel, thirtytwo=False, mom="000"):
    logging.debug("{} {} {} {} {}".format(strangeness, isospin, channel, thirtytwo, mom))
    if thirtytwo:
        basedir = os.path.join(expected_levels_path, "32^3_240/mom_{}".format(mom))
    else:
        basedir = os.path.join(expected_levels_path, "24^3_390/mom_{}".format(mom))

    statistics = "bosonic"
    if args.hadrontype is "baryon":
        statistics = "fermionic"

    expected_isomap = {"0": "I=0", "1": "I=1", "1h": "2I=1", "2": "I=2", "3h": "2h=3", "5h": "2h=5"}
    filename = "{}_{}_S={}_levels.txt".format(statistics, expected_isomap[isospin], strangeness)

    filepath = os.path.join(basedir, filename)
    logging.info("opening {}".format(filepath))
    expected_level_file = open(filepath, "r")
    chan = channel.split("_")[0]
    start = False
    expectedleveltxt = ""
    for line in expected_level_file:
        if " " + chan in line:
            start = True
        if start:
            expectedleveltxt += line
            if not line.strip():
                break

    expected_levels = [level.strip() for level in re.findall("\) (.*)", expectedleveltxt)]
    return expected_levels


def get_ops(args, expected_levels):
    level_num = 0
    already_added = []
    threshold = 10000
    logging.debug(repr(expected_levels))
    for level in expected_levels:
        logging.info(level)
        if args.review:
            raw_input("Look OK?")
        level_num += 1
        if args.threshold and level_num > threshold+3:
            logging.info("Stopping due to threshold")
            args.outfile.write("# Stopping due to threshold\n")
            break
        args.outfile.write("\n# level {} {} \n".format(level_num, level))
        try:
            opset = irreps.translate_name_to_irrep(level)
        except NotImplementedError, e:
            logging.warn("This level {} is not supported {}, skipping".format(level, e))
            args.outfile.write("# level is not supported: {}, skipping \n".format(e))
            threshold = level_num
            # if threshold > 1000:
            #     print "threshold", level_num
            #     args.outfile.write("# Setting threshold to level {}\n".format(level_num))
            #     threshold = level_num
            continue
        except ValueError, e:
            logging.info("level is a single hadron")
            if args.notallsingles:
                ops = irreps.operator_for_singlehadron(level, "PSQ{}".format(sum(int(i)**2 for i in args.momentum)))
                single_hadron(ops)
            continue
        found_one = False
        print opset
        for op in opset:
            p1, p2 = op
            if p1[3] is None or p2[3] is None:
                args.outfile.write("# {} contains particle which primary operator"
                                   " was never defined! \n\n\n".format(level))
                return
                found_one = True
                continue
            mom1,mom2 = re.findall("PSQ([0-9])", level)
            if swap(p1[0], p2[0]):
                logging.info("swapping")
                p1, p2 = p2, p1
                mom1, mom2 = mom2, mom1
            flavor1 = flavor_type(p1[0])
            flavor2 = flavor_type(p2[0])
            #if threshold > 1000 and ((flavor1 != "pion" or flavor2 != "pion") or "rho" in level):
            if threshold > 1000 and level.count("pi") != 2:
                print "threshold", level_num
                args.outfile.write("# Setting threshold to level {}\n".format(level_num))
                threshold = level_num
            if threshold < 1000 and level.count("pi") == 2:
                args.outfile.write("# Found pi pi above threshold. stopping\n".format(level_num))
                threshold = -10

            if flavor2 == "kaon":
                flavor2 = "kbar"
            if (flavor1,p1[2]) == ("eta","A1gp") or (flavor2,p2[2]) == ("eta","A1gp"):
                print "OZI SUPRESSED"
                args.outfile.write("# {} SHOULD BE OZI SUPRESSED!! \n".format(level))
                found_one = True
                continue
            oplines = custom_opline(level, mom1, mom2, p1, p2, flavor1, flavor2, args.channel, args.outfile)
            print oplines
            if oplines is None:
                logging.warn("failed to make this op {} {}".format(p1, p2))
                continue
            logging.debug("success on making op {} {}".format(p1, p2))
            found_one = True
            for opline in oplines:
                if opline in already_added:
                    logging.warn("This operator already added!")
                    args.outfile.write("# This operator already added! \n")
                    args.outfile.write("# {}".format(opline))
                else:
                    args.outfile.write(opline)
                    already_added.append(opline)
                if "eta" in opline:
                    args.outfile.write("# put in ss operators for every uu operator\n")
                    philine = opline.replace("eta", "phi")
                    if philine in already_added:
                        logging.warn("This operator already added!")
                        args.outfile.write("# This operator already added! \n")
                        args.outfile.write("# {}".format(philine))
                    else:
                        args.outfile.write(philine)
                        already_added.append(philine)
        if not found_one:
            logging.critical("Did not find any for this level {} ABORT!!".format(level))
            print flavor1, flavor2
            if (flavor1, flavor2) == ("eta", "pion"):
                logging.critical("This is an eta_pion state we have no ops for!!")
                args.outfile.write("# This is an eta_pion state we have no coeffs for!!\n".format(level_num))
                raw_input("WTF")
            else:
                exit()
    return already_added


def secondary(ops, count):
    args.outfile.write("\n# SECONDARY OPERATORS----------------\n\n")
    newops = []
    for op in ops[:count]:
        words = op.split()
        base, chan, mom1, irrep1, op1, mom2, irrep2, op2 = words
        op1, op2 = op1.strip(']")'), op2.strip(']")')
        _, type1, type2 = base.split("_")
        psqr_1 = sum([int(i)**2 for i in mom1[4:-1].split(",")])
        psqr_2 = sum([int(i)**2 for i in mom2[4:-1].split(",")])
        m1 = irreps.momentums[psqr_1]
        m2 = irreps.momentums[psqr_2]
        isomap = {"pion": "1", "kaon": "1h", "eta": "0", "phi": "0"}
        op1_choices = getprunedops(isomap[type1], m1)[irrep1]
        op2_choices = getprunedops(isomap[type2], m2)[irrep2]
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
    parser.add_argument("-I", "--isospin", choices=["0", "1", "1h", "2"],
                        help="select isospin")
    parser.add_argument("-S", "--strangeness", choices=["0", "1", "2"],
                        help="select strangeness")
    parser.add_argument("-C", "--channel", type=str,
                        help="select channel e.g. A1up_1")
    parser.add_argument("-m", "--momentum", choices=["000", "001", "002", "011", "111"],
                        default="000", help="momentum")
    parser.add_argument("-N", "--number", type=int, required=False,
                        help="number of expected levels to make")
    parser.add_argument("-t", "--threshold", action="store_true", required=False,
                        help="stop after threshold")
    parser.add_argument("--baryon", action="store_true", required=False,
                        help="make baryons instead of mesons")
    parser.add_argument('outstub', nargs='?', type=str,
                        default=None)
    parser.add_argument("-32", "--thirtytwo", action="store_true",
                        default=False, help="use 32^3 expected levels")
    parser.add_argument("--notallsingles", action="store_true",
                        default=False, help="don't include all single hadrons, just the ones in the elevels")
    parser.add_argument("--secondary", type=int, default=0,
                        help="Prompt user for seconday operators,"
                        " the number is how many operators to make secondaries for")

    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(format='# %(levelname)s: %(message)s', level=logging.DEBUG)
        logging.debug("Verbose debuging mode activated")
    else:
        logging.basicConfig(format='# %(levelname)s: %(message)s', level=logging.WARN)

    args.momray = args.momentum.replace("1","+").replace("2","#")
    print args.momray

    get_unspecified_parameters(args)

    psqr = sum(int(i)**2 for i in args.momentum)
    print all_permutations(psqr, nobackies=True)
    for p in all_permutations(psqr, nobackies=True):
        if args.outstub is None:
            args.outfile = sys.stdout
        else:
            args.outfile = open("{}_{}.txt".format(args.outstub, p), 'w')
        args.momray = p.replace("-1", "-").replace("-2","=").replace("1","+").replace("2","#")
        if not args.notallsingles:
            all_single_hadrons(args.momray)

        args.outfile.write("# using mom={}".format(p))

        get_unspecified_parameters(args)
        expected_levels = read_expected_levels(args.strangeness, args.isospin, args.channel, args.thirtytwo, mom=args.momentum)
        logging.info(repr(expected_levels))
        ops = get_ops(args, expected_levels[:args.number])
        if args.secondary:
            secondary(ops, args.secondary)
