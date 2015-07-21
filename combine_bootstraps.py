#!/usr/bin/env python2
import argparse
import logging
import plot_files
import numpy as np
import pandas as pd
import pandas_reader

def writeformat(x):
    """Format complex number into paren grouped format"""
    return [int(x[0]) ,"({},{})".format(np.real(x[1]), np.imag(x[1]))]

def writerealformat(x):
    """Format complex number into paren grouped format"""
    return [int(x[0]) , float(x[1])]

def all_same(items):
    return all(x == items[0] for x in items)


def read_bootstraps(files):

    if "decay_constants" in files[0]:
        data = [pd.read_csv(f,comment='#', names=["decay"]) for f in files]

    else:
        data = [pd.read_csv(f,comment='#', names=["bootstrap", "mass", "amp1", "amp2"]) for f in files]

    return data


def mk_mp(data):
    mk = data[0]["mass"]
    mpi = data[1]["mass"]
    result = 2.0*(mk**2)-(mpi**2)
    return result


def weighted(data, weightfile):
    with open(weightfile) as  weights:
       ws = [float(w) for w in weights.readline().split(",")]
    print ws
    if len(data) != len(ws):
        logging.error("weights and data do not match length!")
        raise RuntimeError("weighted averge requires a weight for each file")

    #result = sum(data[i]["mass"]*ws[i] for i in range(len(data)) )
    result = sum(data[i]*ws[i] for i in range(len(data)) )
    if len(result.columns) > 1:
        result["bootstrap"] = list(range(len(result)))
        result = result.set_index("bootstrap")
    return result


def combine_files(files, options):
    data = read_bootstraps(files)

    data = [df.dropna(axis=1,how='all') for df in data]

    assert all_same([df.shape for df in data]), "They are not all the same size!"

    if options.function == "2mk-mp":
        if len(data) != 2:
            raise RuntimeError("Ratio requires exactly 2 correlators")
        if "ud-s" not in files[0] and "ud-ud" not in files[1]:
            raise RuntimeError("requires mk file then mpi file")
        result = mk_mp(data)
    elif options.function == "weights":
        result = weighted(data, options.weights)
    elif options.function == "ratio":
        result = data[0]/data[1]

    if "decay_constants" in files[0]:
        # outputtext = result["decay"].to_csv(None)
        # print outputtext
        txt = [ "{}".format(i) for i in result["decay"]]
        outputtext = "# shifted\n" + "\n".join(txt)
    else:
        try:
            outputtext = result.to_csv(None, index_label="#bootstrap").replace(",", ", ").replace("amp2", "amp2,    ensamble mean: {}, {}, {}".format(*result.mean()))
        except TypeError:
            outputtext = result.to_csv(None, index_label="#bootstrap").replace(",", ", ")
        except IndexError:
            outputtext = result.to_csv(None, index_label="#bootstrap").replace(",", ", ")


    #outputtext = formated.to_csv(None, sep=",", index=False, header=False, index_label="#time").replace(",", ", ")
    if options.output_stub:
        outfile = "{}".format(options.output_stub)
        with open(outfile, 'w') as ofile:
            ofile.write(outputtext)
            logging.info("wrote ave cor to {}".format(outfile))
    else:
        print outputtext

if __name__ == "__main__":
    functs = ["ratio", "2mk-mp", "weights"]
    parser = argparse.ArgumentParser(description="average data files")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="increase output verbosity")
    parser.add_argument("-o", "--output_stub", type=str, required=False,
                        help="stub of name to write output to")
    parser.add_argument("-f", "--function", type=str, required=True, choices=functs,
                        help="function to use to combine the data")
    parser.add_argument("-n", "--negate", action="store_true",
                        help="write the negative of the correlator")
    parser.add_argument("-w", "--weights", required=False, type=str,
                        help="file with the weights in them")
    parser.add_argument('files', metavar='f', type=str, nargs='+',
                        help='files to plot')
    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.DEBUG)
        logging.debug("Verbose debuging mode activated")
    else:
        logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.INFO)

    if args.negate:
        logging.info("WRITING THE NEGATIVE")

    logging.info(args.files)

    combine_files(args.files, args)
