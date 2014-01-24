#!/usr/bin/env python
import argparse
import logging
import readinput
import os
import re
import irreps
import sys

expected_levels_path = "/home/colin/research/notes/hadron_spectrum/expectedlevels/final_results"
#operators_path = "/latticeQCD/raid6/yjhang/multi_hadron_pruning_operators"
operators_path = "/latticeQCD/raid6/bfahy/operators"

def not_implemented(description, default=""):
    logging.error("No primary operator set for particle {}!".format(description))
    logging.error("Unable to handle particles which have no primary operator set when run in script mode")
    return None

readinput.askoperator = not_implemented



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
        args.isospin = readinput.selectchoices(["0","1","1h"], default="1")

    if not args.strangeness:
        print("Select strangeness")
        args.strangeness = readinput.selectchoices(["0","1","2"], default="0")

    channeldir = os.path.join(operators_path, "BI{I}S{S}".format(I=args.isospin, S=args.strangeness))
    logging.debug("channel dir is {}".format(channeldir))
    channel_list = os.listdir(channeldir)
    if args.channel:
        if args.channel not in channel_list:
            logging.critical("format of input channel is not correct! use e.g. {}".format(channel_list[0]))
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
    for level in expected_levels:
        args.outfile.write("# level {} {} \n".format(level_num, level))
        level_num += 1
        try:
            opset = irreps.translate_name_to_irrep(level)
        except NotImplementedError, e:
            logging.warn("This level {} is not supported {}, skipping".format(level, e))
            args.outfile.write("# level is not supported: {}, skipping \n".format(e))
            continue
        for op in opset:
            p1, p2 = op
            if p1[3] == None or p2[3] == None:
                args.outfile.write("# {} contains particle which primary operator was never defined! \n".format(level))
                continue
            if swap(p1[0], p2[0]):
                logging.info("swapping")
                p1, p2 = p2, p1
            flavor1 = flavor_type(p1[0])
            flavor2 = flavor_type(p2[0])
            if flavor2 == "kaon":
                flavor2 = "kbar"
            flavor_folder = "_".join((flavor1, flavor2))
            opdir = os.path.join(args.opsdir, flavor_folder)
            filename = "S={}_{}_{}_{}_{}_0".format(args.strangeness, p1[1], p1[2], p2[1], p2[2])
            filepath = os.path.join(opdir,filename)
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
                for line in opfile:
                    if re.match(opexpression, line) and re.match(momexpression, line):
                        args.outfile.write(line)
                        if "eta" in line:
                            args.outfile.write("# put in ss operators for every uu operator\n")
                            args.outfile.write(line.replace("eta", "phi"))
            except IOError:
                logging.info("{} didn't exist".format(filename))
                args.outfile.write("# {} didn't exist\n".format(filename))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="compute and plot effective masses")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="increase output verbosity")
    parser.add_argument("-I", "--isospin", choices=["0","1","1h"],
                        help="select isospin")
    parser.add_argument("-S", "--strangeness", choices=["0","1","2"],
                        help="select strangeness")
    parser.add_argument("-C", "--channel", type=str,
                        help="select channel e.g. A1up_1")
    parser.add_argument('outfile', nargs='?', type=argparse.FileType('w'),
                        default=sys.stdout)
    parser.add_argument("-32", "--thirtytwo", action="store_true",
                        default=False, help="use 32^3 expected levels")

    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(format='# %(levelname)s: %(message)s', level=logging.DEBUG)
        logging.debug("Verbose debuging mode activated")
    else:
        logging.basicConfig(format='# %(levelname)s: %(message)s', level=logging.WARN)


    get_unspecified_parameters(args)

    expected_levels = read_expected_levels(args)
    get_ops(args, expected_levels)
