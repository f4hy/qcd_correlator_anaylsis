#!/usr/bin/env python
import logging
import numpy as np
import pandas as pd
from scipy import linalg as LA
from numpy.linalg import cond
import determine_operators
import argparse
import pandas_reader
import os.path

DIAGTOL = 0.009


def compare_matrix(A, B):
    logging.debug("comparing \n{a} to \n{b}".format(a=A, b=B))
    return np.allclose(A, B)


def hermitionize(M):
    """Take a matrix, and return hermition one"""
    return (M+M.H)/2.0


def diagonalize(correlator_pannel, t0, td, generalized=False):
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

    if generalized:
        logging.info("running generalized eigensolver")
        evals, evecs = LA.eigh(A, b=B)  #gerenalized eig problem, eigh works only if hermitian
        evecs = np.matrix(evecs)
        V = evecs
    else:
        # Instead of using generalized eigen problem, we could solve a
        # regular eigen problem involving Binvqrt
        Binvsqrt =  LA.inv(LA.sqrtm(B))
        logging.info("condition number: {}".format(cond(Binvsqrt*A*Binvsqrt)))
        evals, evecs = LA.eigh(Binvsqrt*A*Binvsqrt)
        evecs = np.matrix(evecs)
        V = np.matrix(Binvsqrt)*evecs

    if min(evals) < 0.05:
        logging.warn("Warning, low eigenvalue detected. Eval={}".format(min(evals)))
    else:
        logging.info("lowest eigenvalue={}".format(min(evals)))       
    logging.debug("eigen values are {}".format(evals))
    logging.debug("eigen vectors are {}".format(evecs))
    n = len(evecs)
    logging.debug("Matrix size {N}x{N}".format(N=n))

    def rotate(x):
        M = np.matrix(np.resize(x, (n, n)))
        M = hermitionize(M)
        D = V.H * M * V
        R = np.array(D).flatten()
        return R

    diag = correlator_pannel.apply(rotate, "items")
    diag.items = ["{}{}".format(i,j) for i in reversed(range(n)) for j in reversed(range(n))]

    # This method simultaniously diagaonlizes at t0 and td. Should be
    # identity at t0 and the eigenvalues at td
    assert compare_matrix(np.reshape(diag.major_xs(t0).mean().values, (n, n)),
                          np.identity(n)), "Rotation error: is not ~identity at t0"
    assert compare_matrix(np.reshape(diag.major_xs(td).mean().values, (n, n)),
                          np.diag(evals)), "Rotation error: Is not ~Lambda at td"

    return diag


def parenformat(x):
    """Format complex number into paren grouped format"""
    return "({},{})".format(np.real(x), np.imag(x))


def parenformatn(x):
    """Format list of complex numbers into paren grouped format"""
    return ["({},{})".format(np.real(i), np.imag(i)) for i in x]


def write_cor_matrix(correlator_pannel, outputwild, ops, suffix=""):
    """Write the correlator pannel to files """
    ops = map(str,ops)
    if not os.path.exists(os.path.dirname(outputwild)):
        os.makedirs(os.path.dirname(outputwild))
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
    parser.add_argument("-g", "--generalized", action="store_true",
                        help="used generalized eigen solver")
    parser.add_argument("-a", "--analyize", type=str,
                        help="run effective mass analysis code on results and store in folder")
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
    diag = diagonalize(p, args.tnaught, args.tstar, generalized=args.generalized)
    diagave = diag.mean(2)

    levels = range(len(args.operators))
    if args.outputformat:
        write_cor_matrix(diag, args.outputformat, levels, suffix=".full")
        write_cor_matrix(diagave, args.outputformat, levels, suffix=".ave")

    if args.analyize:
        exe_folder, _ = os.path.split(os.path.realpath(__file__))
        output_dir, outputformat = os.path.split(args.outputformat)
        outputformat += ".full"
        opstring = "-r " + " -r ".join([str(l) for l in levels])
        runstring = '{}/main.py -i {}/ -f {} -nv -o {}/ {}'.format(exe_folder, output_dir,
                                                                   outputformat, args.analyize, opstring)
        logging.info("executing {}".format(runstring))
        os.system(runstring)
