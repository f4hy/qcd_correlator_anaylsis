#!/usr/bin/env python
import correlator
import struct
import build_corr
import logging
import numpy as np
import pandas as pd
from scipy import linalg as LA
import determine_operators
import argparse
import pandas_reader

DIAGTOL = 0.001

def compare_matrix(A, B):
    logging.debug("comparing \n{a} to \n{b}".format(a=A, b=B))
    normed_diff = min(LA.norm(A-B) , LA.norm(np.absolute(A)-B))
    logging.debug("difference of {}".format(normed_diff))
    return normed_diff

def diagonalize(correlator_pannel, t0, ts):
    length = p.shape[0]
    n = int(np.sqrt(length))
    A = np.matrix(np.reshape(p.major_xs(3).values,(n,n)))
    B = np.matrix(np.reshape(p.major_xs(5).values,(n,n)))

    evals, evecs = LA.eigh(A, b=B)
    evecs = np.matrix(evecs)

    def rotate(x):
        n = int(np.sqrt(len(x)))
        logging.debug("Matrix size {N}x{N}".format(N=n))
        M = np.matrix(np.resize(x,(n,n)))
        D = evecs.H * M * evecs
        R = np.array(D).flatten()
        return R

    diag = p.apply(rotate, "items")

    assert (compare_matrix(np.reshape(diag.major_xs(ts).values,(n,n)), np.identity(n)) < DIAGTOL), "Rotation error: is not ~identity at t-star"
    assert (compare_matrix(np.reshape(diag.major_xs(t0).values,(n,n)), np.diag(evals)) < DIAGTOL), "Rotation error: Is not ~Lambda at t0"

    return diag


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Diagonalize correlator matrix")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="increase output verbosity")
    parser.add_argument("-f", "--filewild", type=str, required=True,
                        help="fromat of the correlator files in the directory\n\n"
                        "e.g. cor-snk{}_src{}.dat where {} are replaced with operator strings")
    parser.add_argument("-r", "--operators", action='append', required=False,
                        help="operator to make e.g. -r etap000DDL7Egp1")
    parser.add_argument("-i", "--input-dir", type=str, required=True,
                        help="directory to read files from")
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
            print "Error: no operators found"
            parser.print_help()
            parser.exit()
        args.operators = ops


    for snk in ops:
        for src in ops:
            filename = args.input_dir + args.filewild.format(snk,src)
            logging.info("reading {}".format(filename))
            df = pandas_reader.read_configcols_paraenformat(filename)

            exit()
            pandas_reader.read_full_paraenformat_mean(filename)

    exit()


    logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.DEBUG)
    d00 = {1000: {3:  1.0, 4: 2.0}, 1001: {3: 3.0, 4: 4.0} , 1002: {3: 5.0, 4: 0.0}}
    d01 = {1000: {3: -1.0, 4: 3.0}, 1001: {3: 3.0, 4: 4.0} , 1002: {3: 5.0, 4: 0.0}}
    d10 = {1000: {3:  1.0, 4: 2.0}, 1001: {3: 3.0, 4: 4.0} , 1002: {3: 5.0, 4: 0.0}}
    d11 = {1000: {3: -1.0, 4: 3.0}, 1001: {3: 3.0, 4: 4.0} , 1002: {3: 5.0, 4: 0.0}}

    d00 = {1000: {3:  1.0, 4: 2.0}}
    d01 = {1000: {3: -1.0, 4: 3.0}}
    d10 = {1000: {3:  1.0, 4: 2.0}}
    d11 = {1000: {3: -1.0, 4: 3.0}}

    # moving pion data 0=ss0 1=ss1
    d00 = {1000: {3:  np.complex(25.5398393366745,0.0119055782214339), 5: np.complex(15.0325100472387,0.0158512178206725)}}
    d01 = {1000: {3:  np.complex(-12.2271307837365,-0.0103628313697365), 5: np.complex(-10.20183343045,-0.00381791576319435)}}
    d10 = {1000: {3:  np.complex(-12.2625047443711,0.039971501444209), 5: np.complex(-10.2155376553752,0.0406523146607429)}}
    d11 = {1000: {3:  np.complex(59.1440618668479,0.037089030400901), 5: np.complex(43.9093726008716,0.046347681872425)}}

    # fake moving pion data 0=ss0 1=ss1 forced herm
    d00 = {1000: {3:  np.complex(25.5398393366745,0.0119055782214339), 5: np.complex(15.0325100472387,0.0158512178206725)}}
    d01 = {1000: {3:  np.complex(-12.2271307837365,-0.0103628313697365), 5: np.complex(-10.20183343045,-0.00381791576319435)}}
    d10 = {1000: {3:  np.complex(-12.2271307837365,0.0103628313697365), 5: np.complex(-10.20183343045,0.00381791576319435)}}
    d11 = {1000: {3:  np.complex(59.1440618668479,0.037089030400901), 5: np.complex(43.9093726008716,0.046347681872425)}}


    df00 = pd.DataFrame(d00)
    df01 = pd.DataFrame(d01)
    df10 = pd.DataFrame(d10)
    df11 = pd.DataFrame(d11)

    p = pd.Panel({"00": df00, "01": df01, "10": df10, "11": df11})

    diag = diagonalize(p, 3, 5)
    print np.matrix(np.reshape(diag.major_xs(3).values,(2,2)))
    exit()
