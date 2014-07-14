#!/usr/bin/env python
import logging
import argparse
from collections import namedtuple
import particle_operators
from subductions import subductions

properties = namedtuple('particleproperties', ['I', 'P', 'G', "spin"])

names = {}
# mesonnames
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
names["K0star"]   = properties("1/2", 1, None, 0)
# Baryon names ['I', 'P', 'G', "spin"]
names_baryons = {}
# Baryons have different properties for excited versions
names_baryons["N"]   = properties("1/2", 1, None, "1/2")
names_baryons["N-1440"]   = properties("1/2", 1, None, "1/2")
names_baryons["N-1520"]   = properties("1/2", -1, None, "3/2")
names_baryons["N-1535"]   = properties("1/2", -1, None, "1/2")
names_baryons["N-1650"]   = properties("1/2", -1, None, "1/2")
names_baryons["N-1675"]   = properties("1/2", -1, None, "5/2")
names_baryons["N-1680"]   = properties("1/2", 1, None, "5/2")
names_baryons["N-1700"]   = properties("1/2", -1, None, "3/2")
names_baryons["N-1710"]   = properties("1/2", 1, None, "1/2")
names_baryons["N-1720"]   = properties("1/2", 1, None, "3/2")
names_baryons["N-1860"]   = properites("1/2", 1, None, "5/2")
names_baryons["N-1875"]   = properties("1/2", -1, None, "3/2")
names_baryons["N-1880"]   = properties("1/2", 1, None, "1/2")
names_baryons["N-1895"]   = properties("1/2", -1, None, "1/2")
names_baryons["N-1900"]   = properties("1/2", 1, None, "3/2")
names_baryons["N-1990"]   = properties("1/2", 1, None, "7/2")
names_baryons["N-2000"]   = properties("1/2", 1, None, "5/2")
names_baryons["N-2040"]   = properties("1/2", 1, None, "3/2")
names_baryons["N-2060"]   = properties("1/2", -1, None, "5/2")
names_baryons["N-2100"]   = properties("1/2", 1, None, "1/2")
names_baryons["N-2120"]   = properties("1/2", -1, None, "3/2")
names_baryons["N-2190"]   = properties("1/2", -1, None, "7/2")
names_baryons["N-2220"]   = properties("1/2", 1, None, "9/2")
names_baryons["N-2250"]   = properties("1/2", -1, None, "9/2")
names_baryons["N-2300"]   = properties("1/2", 1, None, "1/2")
names_baryons["N-2570"]   = properties("1/2", -1, None, "5/2")
names_baryons["N-2600"]   = properties("1/2", -1, None, "11/2")
names_baryons["N-2700"]   = properties("1/2", 1, None, "13/2")

# properties of lambdas
names_baryons["Lambda"]   = properties(0, 1, None, "1/2")
names_baryons["Lambda-1405"] = properties(0, -1, None, "1/2")
names_baryons["Lambda-1520"] = properties(0, -1, None, "3/2")
names_baryons["Lambda-1600"] = properties(0, 1, None, "1/2")
names_baryons["Lambda-1670"] = properties(0, -1, None, "1/2")
names_baryons["Lambda-1690"] = properties(0, -1, None, "3/2")
names_baryons["Lambda-1800"] = properties(0, -1, None, "1/2")
names_baryons["Lambda-1810"] = properties(0, 1, None, "1/2")
names_baryons["Lambda-1820"] = properties(0, 1, None, "5/2")
names_baryons["Lambda-1830"] = properties(0, -1, None, "5/2")
names_baryons["Lambda-1890"] = properties(0, 1, None, "3/2")
names_baryons["Lambda-2020"] = properties(0, 1, None, "7/2")
names_baryons["Lambda-2100"] = properties(0, -1, None, "7/2")
names_baryons["Lambda-2110"] = properties(0, 1, None, "5/2")
names_baryons["Lambda-2325"] = properties(0, -1, None, "3/2")
names_baryons["Lambda-2350"] = properties(0, 1, None, "9/2")

# properties of deltas
names_baryons["Delta-1232"] = properties("3/2", 1, None, "3/2")
names_baryons["Delta-1600"] = properties("3/2", 1, None, "3/2")
names_baryons["Delta-1620"] = properties("3/2", -1, None, "1/2")
names_baryons["Delta-1700"] = properties("3/2", -1, None, "3/2")
names_baryons["Delta-1750"] = properties("3/2", 1, None, "1/2")
names_baryons["Delta-1900"] = properties("3/2", -1, None, "1/2")
names_baryons["Delta-1905"] = properties("3/2", 1, None, "5/2")
names_baryons["Delta-1910"] = properties("3/2", 1, None, "1/2")
names_baryons["Delta-1920"] = properties("3/2", 1, None, "3/2")
names_baryons["Delta-1930"] = properties("3/2", -1, None, "5/2")
names_baryons["Delta-1940"] = properties("3/2", -1, None, "3/2")
names_baryons["Delta-1950"] = properties("3/2", 1, None, "7/2")
names_baryons["Delta-2000"] = properties("3/2", 1, None, "5/2")
names_baryons["Delta-2150"] = properties("3/2", -1, None, "1/2")
names_baryons["Delta-2200"] = properties("3/2", -1, None, "7/2")
names_baryons["Delta-2300"] = properties("3/2", 1, None, "9/2")
names_baryons["Delta-2350"] = properties("3/2", -1, None, "5/2")
names_baryons["Delta-2390"] = properties("3/2", 1, None, "7/2")
names_baryons["Delta-2400"] = properties("3/2", -1, None, "9/2")
names_baryons["Delta-2420"] = properties("3/2", 1, None, "11/2")
names_baryons["Delta-2750"] = properties("3/2", -1, None, "13/2")
names_baryons["Delta-2950"] = properties("3/2", 1, None, "15/2")

# properties of xis
names_baryons["Xi-1530"] = properties("1/2", 1, None, "3/2")
names_baryons["Xi-1820"] = properties("1/2", -1, None, "3/2")

# properties of omega
names_baryons["Omega"] = properties(0, 1, None, "3/2")

# properties of sigma
names_baryons["Sigma"] = properties(1, 1, None, "1/2")
names_baryons["Sigma-1385"] = properties(1, 1, None, "3/2")
names_baryons["Sigma-1580"] = properties(1, -1, None, "3/2")
names_baryons["Sigma-1620"] = properties(1, -1, None, "1/2")
names_baryons["Sigma-1660"] = properties(1, 1, None, "1/2")
names_baryons["Sigma-1670"] = properties(1, -1, None, "3/2")
names_baryons["Sigma-1750"] = properties(1, -1, None, "1/2")
names_baryons["Sigma-1770"] = properties(1, 1, None, "1/2")
names_baryons["Sigma-1775"] = properties(1, -1, None, "5/2")
names_baryons["Sigma-1840"] = properties(1, 1, None, "3/2")
names_baryons["Sigma-1880"] = properties(1, 1, None, "1/2")
names_baryons["Sigma-1915"] = properties(1, 1, None, "5/2")
names_baryons["Sigma-1940"] = properties(1, -1, None, "3/2")
names_baryons["Sigma-2000"] = properties(1, -1, None, "1/2")
names_baryons["Sigma-2030"] = properties(1, 1, None, "7/2")
names_baryons["Sigma-2070"] = properties(1, 1, None, "5/2")
names_baryons["Sigma-2080"] = properties(1, 1, None, "3/2")
names_baryons["Sigma-2100"] = properties(1, -1, None, "7/2")


meson_reps = {"A1": (0, 4), "A2": (3, 6), "E": (2, 4, 5, 6), "T1": (1, 3, 4), "T2": (2, 3, 4, 5)}
reps_meson = {0: ["A1"], 4: ["A1"], 1: ["T1"], 2: ["E", "T2"], 3: ["T1", "T2", "A2"],
              4: ["A1", "E", "T1", "T2"]}
baryon_reps = {"G1": ("1/2", "7/2"), "G2": ("5/2", "7/2"), "H": ("3/2", "5/2", "7/2")}
reps_baryon = {"1/2": ["G1"], "3/2": ["H"], "5/2": ["G2", "H"], "7/2": ["G1", "G2", "H"]}

parity = {-1: "u", 1: "g"}
gparity = {-1: "m", 1: "p"}
momentums = {0: "AR", 1: "OA", 2: "PD", 3: "CD", 4: "OA", 5: "PSQ5", 6: "PSQ6", 7: "PSQ7", 8: "PSQ8"}



def irrep_rest_particle(p):
    if ismeson(p):
        logging.info("particle is a meson, gchecking gparity")
        if p.G:
            return ["{}{}{}".format(rep, parity[p.P], gparity[p.G]) for rep in reps_meson[p.spin]]
        else:
            return ["{}{}".format(rep, parity[p.P]) for rep in reps_meson[p.spin]]
    else:
        logging.info("particle is a baryon")
        return ["{}{}".format(rep, parity[p.P]) for rep in reps_baryon[p.spin]]


def getmomint(s):
    return int(s.strip("PSQAB"))

def isbaryon(p):
    return p.spin in reps_baryon.keys()

def ismeson(p):
    return p.spin in reps_meson.keys()


def irrep_moving_particle(p, momentum):
    mom = getmomint(momentum)
    if mom < 1:
        return irrep_rest_particle(p)
    subs = []
    if ismeson(p):
        reps = reps_meson[p.spin]
    if isbaryon(p):
        reps = reps_baryon[p.spin]
    for rep in reps:
        subs.extend(subductions[(rep, parity[p.P], mom)])
    if p.G:
        return ["{}{}".format(movingrep, gparity[p.G]) for movingrep in subs]
    else:
        return ["{}".format(movingrep) for movingrep in subs]


def particle_name(name):
    logging.debug("identifying {}".format(name))
    basename = name.split("-")[0]
    if basename in names.keys():
        return names[basename]
    elif basename in names_baryons.keys():
        return names_baryons[name]
    else:
        logging.critical("I don't know understand the particle {}".format(basename))
        raise NotImplementedError("Particle name unknown")


def operator_for_singlehadron(name, psqr):
    logging.info("Translating {}".format(name))
    operators = particle_operators.particleDatabase()
    p1 = particle_name(name)
    try:
        irreps = irrep_moving_particle(p1, psqr)
    except KeyError:
        logging.critical("I don't know how to do subductions for this momenta")
        raise NotImplementedError("Unsupported momenta")
    ops = {}
    for irrep in irreps:
        op = operators.read_op(name, irrep, getmomint(psqr))
        logging.info("particle1 %s in %s, primary operator is %s", name, irrep, op)
        ops[irrep] = op
    return ops


def translate_name_to_irrep(name):
    logging.info("Translating {}".format(name))
    name = name.replace("KB", "K")
    operators = particle_operators.particleDatabase()
    if "_" not in name:
        logging.info("Is a single hadron?")
        p1 = particle_name(name)
        logging.info("particle %s at rest irreps: %s", name, " ".join(irrep_rest_particle(p1)))
        raise ValueError("is a single hadron")

    if name.count("_") > 5:
        logging.info("Is a three or four hadron state?")
        raise NotImplementedError("3 or more hadrons")

    particle1, momentum1, particle2, momentum2, _, _ = name.split("_")
    mom1 = getmomint(momentum1)
    mom2 = getmomint(momentum2)
    try:
        p1 = particle_name(particle1)
        p2 = particle_name(particle2)
    except KeyError:
        logging.critical("I don't know understand one of the particle names")
        raise NotImplementedError("Particle name unknown")
    logging.info("particle1 %s at rest irreps: %s", particle1, " ".join(irrep_rest_particle(p1)))
    logging.info("particle2 %s at rest irreps: %s", particle2, " ".join(irrep_rest_particle(p2)))

    iso = [i for i in name.split("_") if i.startswith("iso")][0]

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
