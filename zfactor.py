#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
import logging
import argparse
import plot_files
import pandas as pd
import level_identifier
import pandas_reader
import re

OUTPUT = 25
logging.addLevelName(OUTPUT, "OUTPUT")

pair = re.compile(r'\(([^,\)]+),([^,\)]+)\)')

def parse_pair(s):
    return complex(*map(float, pair.match(s).groups()))


def read_coeffs_file(filename):
    txt = plot_files.lines_without_comments(filename)
    parse_pair = plot_files.parse_pair
    df = pd.read_csv(txt, delimiter=' ', names=["id", "identities", "error"],
                     converters={1: parse_pair, 2: parse_pair}, index_col=0)
    return df


def build_cor_mat(corwild, ops, to):
    N = len(ops)
    cormat = np.matrix(np.zeros((N,N)))
    for i,src in enumerate(ops):
        for j,snk in enumerate(ops):
            logging.info("Reading snk:{}, src:{}".format(snk,src))
            raw_c = plot_files.read_file(corwild.format(snk,src))
            df = raw_c
            Cij = df.ix[df['time'] == to, 'correlator']
            cormat[i,j] = np.array(Cij)[0]
    return cormat

def check_ident(v, cormat):
    # check, should be identity
    should_be_identity = (v.H).dot(cormat).dot(v)
    error = np.max(np.abs((should_be_identity - np.identity(len(v)))))
    if error > 1e-4:
        logging.error(u"v^\u2020 C v is different from identity by max {}".format(error))
    else:
        logging.info(u"v^\u2020 C v is different from identity by max {}".format(error))


def read_emasses(filewild, N, t):
    emasses = np.empty(N)
    for level in range(N):
        df = plot_files.read_file(filewild.format(level))
        emasses[level] = np.array(np.real(df.ix[df["time"] == t, "correlator"]))
    return emasses

def normalize_Zs(Zs):
    A = np.array(Zs.values())
    maximums = np.array([max(np.abs(A[:,i])) for i in range(len(Zs[0]))])
    normed = {k: np.abs(values)/maximums for k,values in Zs.iteritems()}
    return normed

def compute_zfactor(corwild, rotfile, emasswild, ops, t0, t, outputstub, maxlevels):
    raw_v = read_coeffs_file(rotfile)
    N = len(ops)
    v = np.matrix(raw_v.identities.values.reshape((N,N))).T
    cormat = build_cor_mat(corwild, ops, t0)
    check_ident(v, cormat)

    levels_to_make = range(min(len(ops),maxlevels))

    emasses = read_emasses(emasswild, len(ops), t)
    Zs = {}
    for level in levels_to_make:
        v_n = v[:,level]
        Zs[level] = [(cormat[j]*(v_n)*np.exp(emasses[level] * t0 * 0.5)).flat[0] for j in range(len(ops)) ]
        print len(Zs[level])
    for level in levels_to_make:
        logging.info("Z_j for level{}: {}".format(level,repr(Zs[level])))

    print Zs.keys()
    print len(Zs[0])
    normalized_Zs = normalize_Zs(Zs)
    for level in levels_to_make:
        logging.info("normed Z_j for level{}: {}".format(level,str(normalized_Zs[level])))

    if(outputstub):
        with open(outputstub+".out", 'w') as outfile:
            outfile.write("# normalized Zfactors\n")
            for level in levels_to_make:
                for j in range(len(ops)):
                    outfile.write("{:d}{:03d} {}\n".format(j+1, level+1, normalized_Zs[level][j]))

    check_sum(Zs, emasses, t)

def check_sum(Zs, emasses, t):
    logging.info("Checking sum")
    for i in range(len(Zs[0])):
        for j in range(len(Zs[0])):
            C = sum((Zs[level][i]*Zs[level][j])*np.exp(-1.0*emasses[level]*t) for level in Zs.keys())
            print "C_{}{}({}) = {}?".format(i, j, t, C)

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Compute the Zfactors from the correlator and diagonalized coeffs")
    parser.add_argument("-ic", "--inputcorrelatorformat", type=str, required=True,
                        help="Correlator file to read from")
    parser.add_argument("-ir", "--inputrotationcoeffs", type=str, required=True,
                        help="rotationcoeffs file to read from")
    parser.add_argument("-ie", "--inputemass", type=str, required=True,
                        help="emass wildcard to read from")
    parser.add_argument("-t", "--time", type=int, required=True,
                        help="time diagonalization was done")
    parser.add_argument("-to", "--tnaught", type=int, required=True,
                        help="t naught, reference time")
    parser.add_argument("-ops", "--operators", type=str, nargs="+", required=True,
                        help="operator strings, order matters!")
    parser.add_argument("-o", "--output_stub", type=str, required=False,
                        help="stub of name to write output to")
    parser.add_argument("-n", "--number", type=int, required=False,
                        help="restrict to a number of levels", default=1000)
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="increase output verbosity")

    args = parser.parse_args()
    if args.verbose:
        logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.DEBUG)
        logging.debug("Verbose debuging mode activated")
    else:
        logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.INFO)

    compute_zfactor(args.inputcorrelatorformat, args.inputrotationcoeffs, args.inputemass, args.operators, args.tnaught, args.time, args.output_stub, args.number)
