#!/usr/bin/env python

"""This contains methods to write out a correlator in binary format
which is expected by colin's analyis code.

"""
import struct
import build_corr
import logging
import math
from itertools import combinations


def write_binary_correlator_matrix(filename, cor_matrix):

    assert all(c1.compatible(c2) for c1, c2 in combinations(cor_matrix.values(), 2))

    size = int(math.sqrt(len(cor_matrix)))

    first = cor_matrix[(1, 1)]

    f = open(filename, 'wb')
    NSTRINGS = 20
    f.write(struct.pack('i', NSTRINGS))  # Nstrings
    f.write("                                    \0")
    f.write("                                    \0")
    f.write("                                    \0")
    f.write("                           beta     \0")  # beta at 27?
    f.write("                           %        \0")  # a_{s}/a_{t}=?
    f.write("                                    \0")
    f.write("                                    \0")
    f.write("                                    \0")
    f.write("                                    \0")
    f.write("                                    \0")
    f.write("                                    \0")
    f.write("                                    \0")
    f.write("                                    \0")
    f.write("               16              36   \0")  # Lattice
    f.write("                                    \0")
    f.write("                                    \0")
    f.write("                                               %s\0" % first.numconfigs)
    f.write("                                               1 \0")
    f.write("                                                 \0")
    f.write("                                    \0")
    f.write(struct.pack('i', size))  # NCOR
    f.write(struct.pack('i', first.numtimes - 1))  # MAX_DT
    f.write(struct.pack('i', first.numtimes))  # T_PERIOD
    f.write(struct.pack('i', first.numconfigs))  # NBINS

    for cfg in first.configs:
        for t in first.times:
            for j in range(1, size + 1):
                for i in range(1, j + 1):
                    f.write(struct.pack('d', cor_matrix[(i, j)].get(config=cfg, time=t)))

    for cfg in first.configs:
        for diag in range(1, size + 1):
            f.write(struct.pack('d', cor_matrix[(diag, diag)].vev1[cfg]))


if __name__ == "__main__":
    logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.DEBUG)

    # exampledir = "/home/bfahy/r3/meyer_results_moreops/stream_1/data/"
    # opfile1 = exampledir +  "opvals_a1pp_0_optype0_op1_diag.dat"
    # opfile2 = exampledir +  "opvals_a1pp_0_optype5_op1_diag.dat"
    # opfile3 = exampledir +  "opvals_a1pp_0_optype11_op1_diag.dat"
    #matrix = build_corr.matrix_from_opfiles([opfile1, opfile2, opfile3])

    ops = ["a1pp_0_optype{}_op1".format(op) for op in [0, 1, 3, 5, 10, 11]]

    directory = "/home/bfahy/r3/effectivemasses/meyer_moreops/total/"
    corformat = "binned_500_{}_{}.cor"
    vevformat = "binned_500_{0}_{0}.vev1"
    matrix = build_corr.matrix_from_cor_and_vev(directory, corformat, vevformat, ops)

    write_binary_correlator_matrix("newop.data1pp", matrix)
