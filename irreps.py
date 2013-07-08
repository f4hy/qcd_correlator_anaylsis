#!/usr/bin/env python
import logging
import argparse
from collections import namedtuple

properties = namedtuple('particleproperties', ['I', 'P', 'G', "spin"])

names = {}
names["pi"]       = properties(1, -1, -1, 0)
names["a0"]       = properties(1, 1, -1, 0)
names["a1"]       = properties(1, 1, -1, 1)
names["rho"]      = properties(1, -1, 1, 1)
names["b1"]       = properties(1, 1, 1, 1)
names["h1"]       = properties(0, 1, -1, 1)
names["omega"]    = properties(0, -1, -1, 1)
names["f1"]       = properties(0, 1, 1, 1)
names["f2"]       = properties(0, 1, 1, 2)
names["eta"]      = properties(0, -1, 1, 0)
names["etaprime"] = properties(0, -1, 1, 0)
names["K"]        = properties("1/2", -1, None, 0)
names["KB"]       = properties("1/2", -1, None, 0)
names["Kstar"]    = properties("1/2", -1, None, 1)
names["K1"]    = properties("1/2", 1, None, 1)

meson_reps = {"A1": (0, 4), "A2": (3, 6), "E": (2, 4, 5, 6), "T1": (1, 3, 4), "T2": (2, 3, 4, 5)}
reps_meson = {0: ["A1"], 4: ["A1"], 1: ["T1"], 2: ["E", "T2"], 3: ["T1", "T2", "A2"],
              4: ["A1", "E", "T1", "T2"]}
baryon_reps = {"G1": ("1/2", "7/2"), "G2": ("5/2", "7/2"), "H": ("3/2", "5/2", "7/2")}

parity = {-1: "u", 1: "g"}
gparity = {-1: "m", 1: "p"}
momentums = {0: "AR", 1: "OA", 2: "PD", 3: "CD"}

subductions = {}
subductions[("A1", "g", 1)] = ["A1"]
subductions[("A1", "u", 1)] = ["A2"]
subductions[("A2", "g", 1)] = ["B1"]
subductions[("A2", "u", 1)] = ["B2"]
subductions[("E", "g", 1)] = ["A1", "B1"]
subductions[("E", "u", 1)] = ["A2", "B2"]
subductions[("T1", "g", 1)] = ["A2", "E"]
subductions[("T1", "u", 1)] = ["A1", "E"]
subductions[("T2", "g", 1)] = ["B2", "E"]
subductions[("T2", "u", 1)] = ["B1", "E"]

subductions[("A1", "g", 2)] = ["A1"]
subductions[("A1", "u", 2)] = ["A2"]
subductions[("A2", "g", 2)] = ["A2"]
subductions[("A2", "u", 2)] = ["A1"]
subductions[("E", "g", 2)] = ["E"]
subductions[("E", "u", 2)] = ["E"]
subductions[("T1", "g", 2)] = ["A2", "E"]
subductions[("T1", "u", 2)] = ["A1", "E"]
subductions[("T2", "g", 2)] = ["A1", "E"]
subductions[("T2", "u", 2)] = ["A2", "E"]

subductions[("A1", "g", 3)] = ["A1"]
subductions[("A1", "u", 3)] = ["A2"]
subductions[("A2", "g", 3)] = ["B2"]
subductions[("A2", "u", 3)] = ["B1"]
subductions[("E", "g", 3)] = ["A1", "B2"]
subductions[("E", "u", 3)] = ["A2", "B1"]
subductions[("T1", "g", 3)] = ["A2", "B1", "B2"]
subductions[("T1", "u", 3)] = ["A1", "B1", "B2"]
subductions[("T2", "g", 3)] = ["A1", "A2", "B1"]
subductions[("T2", "u", 3)] = ["A1", "A2", "B2"]


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
    logging.info("Translateing {}".format(name))

    particle1, momentum1, particle2, momentum2, _, _ = name.split("_")
    mom1 = int(momentum1[-1])
    mom2 = int(momentum2[-1])
    p1 = particle_name(particle1)
    p2 = particle_name(particle2)
    logging.info("particle1 %s atrest irreps: %s", particle1, " ".join(irrep_rest_particle(p1)))
    logging.info("particle2 %s atrest irreps: %s", particle2, " ".join(irrep_rest_particle(p2)))

    if mom1 > 0:
        logging.info("particle1 %s moving irreps: %s", particle1,
                     " ".join(irrep_moving_particle(p1, momentum1)))

    if mom2 > 0:
        logging.info("particle2 %s moving irreps: %s", particle2,
                     " ".join(irrep_moving_particle(p2, momentum2)))

    for i in irrep_moving_particle(p1, momentum1):
        for j in irrep_moving_particle(p2, momentum2):
            logging.info("operatorfile: {}_{}_{}_{}".format(momentums[mom1], i, momentums[mom2], j))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Determine lattice irreps from particle names")
    parser.add_argument('name', metavar='n', type=str, help='files to plot')
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="increase output verbosity")

    args = parser.parse_args()
    if args.verbose:
        logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.DEBUG)
        logging.debug("Verbose debuging mode activated")
    else:
        logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.INFO)

    translate_name_to_irrep(args.name)
