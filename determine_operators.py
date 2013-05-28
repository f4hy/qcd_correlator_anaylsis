#!/usr/bin/env python
import fnmatch
import logging
import os
import re


def matching_operators(directory, pattern):
    modifiled_pattern = pattern.replace("{}", "*")

    regex = fnmatch.translate(modifiled_pattern)
    print regex
    regex_grouped = regex.replace(".*", "(.*)")
    print regex_grouped

    reobj = re.compile(regex_grouped)
    ops = set()
    for f in os.listdir(directory):
        try:
            ops.add(reobj.match(f).group(1))
            ops.add(reobj.match(f).group(2))
        except AttributeError:
            logging.info("file {} did not match pattern".format(f))
    return ops


def matching_files(directory, pattern):
    modifiled_pattern = pattern.replace("{}", "*")
    return fnmatch.filter(os.listdir(directory), modifiled_pattern)


if __name__ == "__main__":
    logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.DEBUG)
    directory = "/home/bfahy/r2/pruning/special/atrestpions/correlators/"
    form = "corsnk-{}_src-{}.dat"

    ops = matching_operators(directory, form)
    print ops
    files = matching_files(directory, form)
    print files
