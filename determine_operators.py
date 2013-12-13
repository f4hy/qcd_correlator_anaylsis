#!/usr/bin/env python
import fnmatch
import argparse
import logging
import os
import re


def matching_operators(directory, pattern):
    modifiled_pattern = pattern.replace("{}", "*")

    regex = fnmatch.translate(modifiled_pattern)
    # print regex
    regex_grouped = regex.replace(".*", "(.*)")
    # print regex_grouped

    reobj = re.compile(regex_grouped)
    ops = set()
    for f in os.listdir(directory):
        try:
            ops.add(reobj.match(f).group(1))
            ops.add(reobj.match(f).group(2))
        except AttributeError:
            logging.info("file {} did not match pattern".format(f))
        except IndexError:
            logging.debug("doesn't have two wildcards")
    return sorted(ops)


def matching_files(directory, pattern):
    modifiled_pattern = pattern.replace("{}", "*")
    return fnmatch.filter(os.listdir(directory), modifiled_pattern)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="compute and plot effective masses")
    parser.add_argument("-f", "--format", type=str, required=False,
                        help="fromat of the correlator files in the directory\n\n"
                        "e.g. {}-{}.A1gp.conn.dat where {} are replaced with operator strings"
                        "Defaults to 'corsnk-{}_src-{}.dat'",
                        default="corsnk-{}_src-{}.dat")
    parser.add_argument("-i", "--input-dir", type=str, required=True,
                        help="directory to read files from")
    logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.DEBUG)

    args = parser.parse_args()

    directory = args.input_dir
    form = "corsnk-{}_src-{}.dat"

    ops = matching_operators(directory, args.format)
    print " ".join(ops)
    # files = matching_files(directory, args.format)
    # print files
