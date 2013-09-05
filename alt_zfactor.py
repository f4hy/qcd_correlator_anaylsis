#!/usr/bin/env python
import numpy as np
import logging
import argparse
import plot_files
import pandas as pd
import zfactor


def build_cor_mat_error(corwild, ops, to):
    N = len(ops)
    cormat = np.matrix(np.zeros((N, N)), dtype=np.complex128)
    errmat = np.matrix(np.zeros((N, N)), dtype=np.complex128)
    for col, src in enumerate(ops):
        for row, snk in enumerate(ops):
            logging.debug("Reading snk:{}, src:{}".format(snk, src))
            raw_c = plot_files.read_file(corwild.format(snk, src))
            df = raw_c
            Cij = df.ix[df['time'] == to, 'correlator']
            cormat[row, col] = np.array(Cij)[0]
            Eij = df.ix[df['time'] == to, 'error']
            errmat[row, col] = np.array(Eij)[0]
    return cormat, errmat


def read_zrots(filename):
    txt = plot_files.lines_without_comments(filename)
    df = pd.read_csv(txt, delimiter=',', names=["level", "amp", "error"], index_col=0)
    return df


def read_level_fits(filename):
    """
    Read in fit values to the levels. Fits should be single exp
    Should be in the format
    # Level, Amp, Error(Amp), Mass, Error(Mass)
    0, 1.74082463996023, 0.0154835440133309, 0.153284504901198, 0.000705009866455281

    final column chi^2 optional
    """

    txt = plot_files.lines_without_comments(filename)
    df = pd.read_csv(txt, delimiter=',', names=["level", "amp", "amp_error", "mass", "mass_error", "chisqr"], index_col=0)
    return df


def alt_zfactor(corwild, zrotfile, rotfile, ops, t0, outputstub,
                maxlevels=None, normalize=False, reconstruct_stub=None, inputemass=None):
    # zrots = read_zrots(zrotfile)
    fit_values = read_level_fits(zrotfile)
    N = len(fit_values)
    N = len(ops)
    levels_to_make = range(min(N, maxlevels))
    raw_v = zfactor.read_coeffs_file(rotfile)
    v = np.matrix(raw_v.identities.values.reshape((N, N))).T
    roterror = np.matrix(raw_v.error.values.reshape((N, N))).T
    cormat, errmat = build_cor_mat_error(corwild, ops, t0)

    Zs = {}
    err = {}
    for level in levels_to_make:
        zr = fit_values.amp[level]
        #err[level] = zrots.error[level]*np.ones(N)
        for op in range(N):
            v_n = (v[:, level])
            ev_n = np.ravel(roterror[:, level])
            Zs[level] = [np.abs((cormat[j]*(v_n)).flat[0])*np.sqrt(zr) for j in range(len(ops))]
            err[level] = [np.abs((cormat[j]*(v_n)).flat[0])*np.sqrt(fit_values.amp_error[level]) +
                          np.abs((errmat[j]*(v_n)).flat[0])*np.sqrt(zr)
                          for j in range(len(ops))]
            #err[level] = np.sqrt(fit_values.amp_error[level])*np.abs(v_n)+np.sqrt(zr)*np.abs(ev_n)
    normalized_Zs = zfactor.normalize_Zs(Zs, normalize)
    A = np.array(Zs.values())
    maximums = np.array([max(np.abs(A[:, i])) for i in range(len(Zs[0]))])
    if normalize:
        normalized_Zs = {k: np.abs(values)/maximums for k, values in Zs.iteritems()}
        normalized_err = {k: np.abs(values)/maximums for k, values in err.iteritems()}
    else:
        normalized_Zs = Zs
        normalized_err = err

    # print err
    # print normalized_err

    if(outputstub):
        with open(outputstub+".out", 'w') as outfile:
            outfile.write("# normalized Zfactors\n")
            for level in levels_to_make:
                for j in range(N):
                    outfile.write("{:d}{:03d} {} {}\n".format(j+1, level+1,
                                                              normalized_Zs[level][j], normalized_err[level][j]))

    if(reconstruct_stub):
        reconstructed_correaltors(Zs, err, fit_values, ops, reconstruct_stub)


def reconstructed_correaltors(Zs, error, fit_values, ops, stub):
    emasses = fit_values.mass
    emasses_err = fit_values.mass_error
    for i in range(len(Zs[0])):
        for j in range(len(Zs[0])):
            with open("{}.{}.{}.cor".format(stub, ops[i], ops[j]), "w") as outfile:
                for t in range(40):
                    C = sum((Zs[level][i]*np.conj(Zs[level][j]))*np.exp(-1.0*emasses[level]*t) for level in Zs.keys())
                    Cerr = sum((error[level][i]*np.conj(Zs[level][j]))*np.exp(-1.0 * emasses[level] * t)+
                               (Zs[level][i]*np.conj(error[level][j]))*np.exp(-1.0 * emasses[level] * t)+
                               (Zs[level][i]*np.conj(Zs[level][j]))*np.exp(-1.0*emasses[level]*t)*(-1.0*emasses_err[level]*t)
                               for level in Zs.keys())
                    outfile.write("{} ({},{}) ({},{})\n".format(t, np.real(C), np.imag(C), np.real(Cerr), np.imag(Cerr)))


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Compute the Zfactors from the correlator and diagonalized coeffs")
    parser.add_argument("-zr", "--z-rot", type=str, required=True,
                        help="file to read the zrot from")
    parser.add_argument("-ic", "--inputcorrelatorformat", type=str, required=True,
                        help="Correlator file to read from")
    parser.add_argument("-ops", "--operators", type=str, nargs="+", required=True,
                        help="operator strings, order matters!")
    parser.add_argument("-ir", "--inputrotationcoeffs", type=str, required=True,
                        help="rotationcoeffs file to read from")
    parser.add_argument("-t0", "--tnaught", type=int, required=True,
                        help="t naught, reference time")
    parser.add_argument("-o", "--output_stub", type=str, required=False,
                        help="stub of name to write output to")
    parser.add_argument("-r", "--reconstruct_stub", type=str, required=False,
                        help="stub for reconstrcuting the correlators")
    parser.add_argument("-n", "--number", type=int, required=False,
                        help="restrict to a number of levels", default=1000)
    parser.add_argument("-norm", "--normalize", action="store_true", required=False,
                        help="normalized the zfactors")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="increase output verbosity")

    args = parser.parse_args()
    if args.verbose:
        logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.DEBUG)
        logging.debug("Verbose debuging mode activated")
    else:
        logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.INFO)

        alt_zfactor(args.inputcorrelatorformat, args.z_rot, args.inputrotationcoeffs, args.operators,
                    args.tnaught, args.output_stub, args.number, args.normalize, args.reconstruct_stub)
