#!/usr/bin/env python
import logging
import argparse
from collections import namedtuple
import particle_operators

properties = namedtuple('particleproperties', ['I', 'P', 'G', "spin"])

names = {}
names["pi"]       = properties(1, -1, -1, 0)
names["pi1"]      = properties(1, -1, -1, 1)
names["pi2"]      = properties(1, -1, -1, 2)
names["a0"]       = properties(1, 1, -1, 0)
names["a1"]       = properties(1, 1, -1, 1)
names["a2"]       = properties(1, 1, -1, 2)
names["rho"]      = properties(1, -1, 1, 1)
names["rho3"]     = properties(1, -1, 1, 3)
names["b1"]       = properties(1, 1, 1, 1)
names["h1"]       = properties(0, 1, -1, 1)
names["omega"]    = properties(0, -1, -1, 1)
names["omega3"]   = properties(0, -1, -1, 3)
names["f0"]       = properties(0, 1, 1, 0)
names["f1"]       = properties(0, 1, 1, 1)
names["f2"]       = properties(0, 1, 1, 2)
names["f2prime"]  = properties(0, 1, 1, 2)
names["eta"]      = properties(0, -1, 1, 0)
names["etaprime"] = properties(0, -1, 1, 0)
names["eta2"]     = properties(0, -1, 1, 2)
names["phi"]      = properties(0, -1, -1, 1)
names["phi3"]     = properties(0, -1, -1, 3)
names["K"]        = properties("1/2", -1, None, 0)
names["KB"]       = properties("1/2", -1, None, 0)
names["Kstar"]    = properties("1/2", -1, None, 1)
names["KBstar"]   = properties("1/2", -1, None, 1)
names["K1"]       = properties("1/2", 1, None, 1)
names["KB1"]      = properties("1/2", 1, None, 1)
names["K2"]       = properties("1/2", -1, None, 2)
names["K2star"]   = properties("1/2", 1, None, 2)
names["K3star"]   = properties("1/2", 1, None, 3)
names["K0star"]   = properties("1/2", 1, None, 2)

meson_reps = {"A1": (0, 4), "A2": (3, 6), "E": (2, 4, 5, 6), "T1": (1, 3, 4), "T2": (2, 3, 4, 5)}
reps_meson = {0: ["A1"], 4: ["A1"], 1: ["T1"], 2: ["E", "T2"], 3: ["T1", "T2", "A2"],
              4: ["A1", "E", "T1", "T2"]}
baryon_reps = {"G1": ("1/2", "7/2"), "G2": ("5/2", "7/2"), "H": ("3/2", "5/2", "7/2")}

parity = {-1: "u", 1: "g"}
gparity = {-1: "m", 1: "p"}
momentums = {0: "AR", 1: "OA", 2: "PD", 3: "CD", 4: "OA", 5: "PSQ5", 6: "PSQ6"}

subductions = {}
subductions[("A1", "g", 1)] = ["A1"]
subductions[("A1", "u", 1)] = ["A2"]
subductions[("A2", "g", 1)] = ["B1"]
subductions[("A2", "u", 1)] = ["B2"]
subductions[("E", "g", 1)]  = ["A1", "B1"]
subductions[("E", "u", 1)]  = ["A2", "B2"]
subductions[("T1", "g", 1)] = ["A2", "E"]
subductions[("T1", "u", 1)] = ["A1", "E"]
subductions[("T2", "g", 1)] = ["B2", "E"]
subductions[("T2", "u", 1)] = ["B1", "E"]

subductions[("A1", "g", 2)] = ["A1"]
subductions[("A1", "u", 2)] = ["A2"]
subductions[("A2", "g", 2)] = ["B2"]
subductions[("A2", "u", 2)] = ["B1"]
subductions[("E", "g", 2)]  = ["A1", "B2"]
subductions[("E", "u", 2)]  = ["A2", "B1"]
subductions[("T1", "g", 2)] = ["A2", "B1", "B2"]
subductions[("T1", "u", 2)] = ["A1", "B1", "B2"]
subductions[("T2", "g", 2)] = ["A1", "A2", "B1"]
subductions[("T2", "u", 2)] = ["A1", "A2", "B2"]

subductions[("A1", "g", 3)] = ["A1"]
subductions[("A1", "u", 3)] = ["A2"]
subductions[("A2", "g", 3)] = ["A2"]
subductions[("A2", "u", 3)] = ["A1"]
subductions[("E", "g", 3)]  = ["E"]
subductions[("E", "u", 3)]  = ["E"]
subductions[("T1", "g", 3)] = ["A2", "E"]
subductions[("T1", "u", 3)] = ["A1", "E"]
subductions[("T2", "g", 3)] = ["A1", "E"]
subductions[("T2", "u", 3)] = ["A2", "E"]

subductions[("A1", "g", 4)] = ["A1"]
subductions[("A1", "u", 4)] = ["A2"]
subductions[("A2", "g", 4)] = ["B1"]
subductions[("A2", "u", 4)] = ["B2"]
subductions[("E", "g", 4)]  = ["A1", "B1"]
subductions[("E", "u", 4)]  = ["A2", "B2"]
subductions[("T1", "g", 4)] = ["A2", "E"]
subductions[("T1", "u", 4)] = ["A1", "E"]
subductions[("T2", "g", 4)] = ["B2", "E"]
subductions[("T2", "u", 4)] = ["B1", "E"]

subductions[("A1", "u", 5)] = ["A2"]
subductions[("A1", "u", 6)] = ["A2"]


def irrep_rest_particle(p):
    if p.G:
        return ["{}{}{}".format(rep, parity[p.P], gparity[p.G]) for rep in reps_meson[p.spin]]
    else:
        return ["{}{}".format(rep, parity[p.P]) for rep in reps_meson[p.spin]]


def irrep_moving_particle(p, momentum):
    mom = int(momentum[-1])
    if mom < 1:
        return irrep_rest_particle(p)
    subs = []
    for rep in reps_meson[p.spin]:
        subs.extend(subductions[(rep, parity[p.P], mom)])
    if p.G:
        return ["{}{}".format(movingrep, gparity[p.G]) for movingrep in subs]
    else:
        return ["{}".format(movingrep) for movingrep in subs]


def particle_name(name):
    logging.debug("identifying {}".format(name))
    return names[name.split("-")[0]]


def translate_name_to_irrep(name):
    logging.info("Translating {}".format(name))
    name = name.replace("KB","K")
    if "_" not in name:
        logging.info("Is a single hadron?")
        p1 = particle_name(name)
        logging.info("particle %s at rest irreps: %s", name, " ".join(irrep_rest_particle(p1)))
        raise NotImplementedError("is a single hadron")

    if name.count("_") > 5:
        logging.info("Is a three or four hadron state?")
        raise NotImplementedError("3 or more hadrons")

    particle1, momentum1, particle2, momentum2, _, _ = name.split("_")
    mom1 = int(momentum1[-1])
    mom2 = int(momentum2[-1])
    try:
        p1 = particle_name(particle1)
        p2 = particle_name(particle2)
    except KeyError:
        logging.critical("I don't know understand one of the particle names")
        raise NotImplementedError("Particle name unknown")
    logging.info("particle1 %s at rest irreps: %s", particle1, " ".join(irrep_rest_particle(p1)))
    logging.info("particle2 %s at rest irreps: %s", particle2, " ".join(irrep_rest_particle(p2)))

    iso = [i for i in name.split("_") if i.startswith("iso")][0]

    operators = particle_operators.particleDatabase()

    try:
        irreps1 = irrep_moving_particle(p1, momentum1)
        irreps2 = irrep_moving_particle(p2, momentum2)
    except KeyError:
        logging.critical("I don't know how to do subductions for this momenta")
        raise NotImplementedError("Unsupported momenta")

    if mom1 > 0:
        logging.info("particle1 %s moving irreps: %s", particle1,
                     " ".join(irreps1))

    if mom2 > 0:
        logging.info("particle2 %s moving irreps: %s", particle2,
                     " ".join(irreps2))

    for i in irreps1:
        for j in irreps2:
            logging.info("operator file: {}_{}_{}_{}".format(momentums[mom1], i, momentums[mom2], j))

    for irrep in irreps1:
        op = operators.read_op(particle1, irrep, mom1)
        logging.info("particle1 %s in %s, primary operator is %s", particle1, irrep, op)

    for irrep in irreps2:
        op = operators.read_op(particle2, irrep, mom2)
        logging.info("particle2 %s in %s, primary operator is %s", particle2, irrep, op)

    opset = []
    for i in irreps1:
        op1 = operators.read_op(particle1, i, mom1)
        for j in irreps2:
            op2 = operators.read_op(particle2, j, mom2)
            opset.append(((p1, momentums[mom1], i, op1), (p2, momentums[mom2], j, op2)))
            opfile = "{}_{}_{}_{}".format(momentums[mom1], i, momentums[mom2], j)
            if op1 == op2:
                logging.info("cat S\=0_{opfile}_* | grep -r {op1}.*{op2}".format(opfile=opfile, op1=op1, op2=op2))
            else:
                logging.info("cat S\=0_{opfile}_* | grep {op1} |grep {op2}".format(opfile=opfile, op1=op1, op2=op2))
            #print "op1:", "p={}".format(mom1), i, op1, "op2:", "p={}".format(mom2), j, op2

    return opset

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Determine lattice irreps from particle names")
    parser.add_argument('name', metavar='expected_level', type=str,
                        help='expected level to determine irreps from')
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="increase output verbosity")

    args = parser.parse_args()
    if args.verbose:
        logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.DEBUG)
        logging.debug("Verbose debugging mode activated")
    else:
        logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.INFO)

    translate_name_to_irrep(args.name)
