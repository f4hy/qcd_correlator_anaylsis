#!/usr/bin/env python
import logging
import numpy as np
import pandas as pd
from scipy import linalg as LA
import determine_operators
import argparse
import pandas_reader
from level_identifier import readops
from diagonalize import hermitionize


def write_eigenvalues(evalues):
    logging.info("Writing eigen values to {}".format(args.output))
    maxeval = max(evalues)
    with open(args.output, "w") as efile:
        efile.write("#eigenvalue, condition number\n")
        for e in evalues:
            efile.write("{}, {}\n".format(e, maxeval/abs(e)))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Diagonalize correlator matrix")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="increase output verbosity")
    parser.add_argument("-g", "--generalized", action="store_true",
                        help="used generalized eigen solver")
    parser.add_argument("-p", "--principle", action="store_true",
                        help="compute princple correlator (diag on each t)")
    parser.add_argument("-a", "--analyize", type=str,
                        help="run effective mass analysis code on results and store in folder")
    parser.add_argument("-f", "--filewild", type=str, required=True,
                        help="fromat of the correlator files in the directory\n\n"
                        "e.g. cor-snk{}_src{}.dat where {} are replaced with operator strings")
    parser.add_argument("-t", "--time", type=int, required=False,
                        help="t naught, reference time")
    parser.add_argument("-r", "--operators", type=str, required=False,
                        help="file with operators to use")
    parser.add_argument("-i", "--input-dir", type=str, required=True,
                        help="directory to read files from")
    parser.add_argument("-e", "--write_eigenvalues", type=str, required=False,
                        help="just write the eigenvalues")
    parser.add_argument("-m", "--returnmaxeigen", action="store_true", required=False,
                        help="return max eigenvalue")
    parser.add_argument("-o", "--output", type=str, required=False,
                        help="file to write output to")
    parser.add_argument("-c", "--conditionnumber", type=float, required=False,
                        help="condition number to determine cutoff")
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
    else:
        ops = readops(args.operators)
        logging.debug("found operators: {}".format(','.join(ops)))
        args.operators = ops

    cor_matrix = {}
    for snk in args.operators:
        for src in args.operators:
            filename = args.input_dir + args.filewild.format(snk, src)
            logging.info("reading {}".format(filename))
            cor_matrix[snk+src] = pandas_reader.read_configcols_paraenformat(filename)

    correlator_pannel = pd.Panel(cor_matrix)

    length = correlator_pannel.shape[0]
    n = int(np.sqrt(length))
    B = np.matrix(np.reshape(correlator_pannel.major_xs(args.time).mean().values, (n, n)))
    B = hermitionize(B)

    evals = LA.eigvalsh(B)
    logging.info("eigenvalues are {}".format(",".join(map(str, evals))))
    if args.output:
        write_eigenvalues(evals)

    meval = max(evals)
    logging.info("max eigen {}".format(meval))

    if args.conditionnumber:
        cutoff = max(evals)/args.conditionnumber
        logging.info("cutoff to ensure condition number <{} is {}".format(args.conditionnumber,
                                                                          cutoff))
        logging.info("would remove {} states".format(len([i for i in evals if i < cutoff])))
        print cutoff
    else:
        logging.info("returning max eigen")
        print meval
