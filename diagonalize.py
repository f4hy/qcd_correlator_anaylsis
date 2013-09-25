#!/usr/bin/env python
import logging
import numpy as np
import pandas as pd
from scipy import linalg as LA
import determine_operators
import argparse
import pandas_reader

DIAGTOL = 0.009


def compare_matrix(A, B):
    logging.debug("comparing \n{a} to \n{b}".format(a=A, b=B))
    return np.allclose(A, B)


def hermitionize(M):
    """Take a matrix, and return hermition one"""
    return (M+M.H)/2.0


def diagonalize(correlator_pannel, t0, td):
    length = correlator_pannel.shape[0]
    n = int(np.sqrt(length))
    # Here we access the pannel major_xs gives time(n), mean incase it
    # was a multi correlator should have no effect on an already averaged one
    A = np.matrix(np.reshape(correlator_pannel.major_xs(td).mean().values, (n, n)))
    B = np.matrix(np.reshape(correlator_pannel.major_xs(t0).mean().values, (n, n)))
    # Require A and B to be hermition for our generalized eigen value
    # problem method to work. Here we force the matricies to be
    # hermtion. This is justified since it is just useing the other
    # measurement of the same value and averaging them.
    A = hermitionize(A)
    B = hermitionize(B)
    logging.debug("A = {} \n B = {}".format(A, B))

    # Instead of using generalized eigen problem, we could solve a
    # regular eigen problem involving Binvqrt
    # Binvsqrt =  LA.inv(LA.sqrtm(B))
    # evals, evecs = LA.eigh(Binvsqrt*A*Binvsqrt)
    # evecs = np.matrix(evecs)
    # V = np.matrix(Binvsqrt)*evecs

    evals, evecs = LA.eigh(A, b=B)  # eigh works only if hermitian
    evecs = np.matrix(evecs)
    logging.debug("eigen values are {}".format(evals))
    logging.debug("eigen vectors are {}".format(evecs))
    n = len(evecs)
    logging.debug("Matrix size {N}x{N}".format(N=n))

    def rotate(x):
        M = np.matrix(np.resize(x, (n, n)))
        #M = hermitionize(M)
        D = evecs.H * M * evecs
        R = np.array(D).flatten()
        return R

    diag = correlator_pannel.apply(rotate, "items")

    assert compare_matrix(np.reshape(diag.major_xs(t0).mean().values, (n, n)),
                          np.identity(n)), "Rotation error: is not ~identity at t0"
    assert compare_matrix(np.reshape(diag.major_xs(td).mean().values, (n, n)),
                          np.diag(evals)), "Rotation error: Is not ~Lambda at td"

    return diag


def parenformat(x):
    return "({},{})".format(np.real(x), np.imag(x))


def parenformatn(x):
    return ["({},{})".format(np.real(i), np.imag(i)) for i in x]


def write_cor_matrix(correlator_pannel, outputwild, ops, suffix=""):
    for snk in ops:
        for src in ops:
            filename = outputwild.format(snk, src)+suffix
            first = True
            cor = correlator_pannel[snk+src]
            if cor.ndim == 1:
                cor.apply(parenformat).to_csv(filename, sep=" ", header=True, index_label="#time")
                continue
            for n in cor:
                if first:
                    cor[n].apply(parenformat).to_csv(filename, sep=" ", header=True, index_label="#time")
                    first = False
                else:
                    cor[n].apply(parenformat).to_csv(filename, sep=" ", header=False, mode="a")
    logging.info("Wrote correlator matrix to {}{}".format(outputwild.format("SNK", "SRC"), suffix))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Diagonalize correlator matrix")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="increase output verbosity")
    parser.add_argument("-f", "--filewild", type=str, required=True,
                        help="fromat of the correlator files in the directory\n\n"
                        "e.g. cor-snk{}_src{}.dat where {} are replaced with operator strings")
    parser.add_argument("-ts", "--tstar", type=int, required=True,
                        help="time diagonalization was done")
    parser.add_argument("-to", "--tnaught", type=int, required=True,
                        help="t naught, reference time")
    parser.add_argument("-r", "--operators", action='append', required=False,
                        help="operator to make e.g. -r etap000DDL7Egp1")
    parser.add_argument("-i", "--input-dir", type=str, required=True,
                        help="directory to read files from")
    parser.add_argument("-o", "--outputformat", type=str, required=False,
                        help="format to write output")
    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.DEBUG)
        logging.debug("Verbose debuging mode activated")
    else:
        logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.INFO)

    if not args.operators:
        logging.info("Operators not specified, attempting to automagically determine")
        ops = determine_operators.matching_operators(args.input_dir, args.filewild)
        logging.info("Operators automagically found to be {}".format(",".join(ops)))
        if not ops:
            logging.error("Error: no operators found")
            parser.print_help()
            parser.exit()
        args.operators = ops

    cor_matrix = {}
    cor_matrix_multi = {}
    cor_matrix_ave = {}
    for snk in args.operators:
        for src in args.operators:
            filename = args.input_dir + args.filewild.format(snk, src)
            logging.info("reading {}".format(filename))
            cor_matrix[snk+src] = pandas_reader.read_configcols_paraenformat(filename)

    p = pd.Panel(cor_matrix)
    diag = diagonalize(p, args.tnaught, args.tstar)
    diagave = diag.mean(2)

    if args.outputformat:
        write_cor_matrix(diag, args.outputformat, args.operators, suffix=".full")
        write_cor_matrix(diagave, args.outputformat, args.operators, suffix=".ave")
