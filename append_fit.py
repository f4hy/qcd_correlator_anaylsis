#!/usr/bin/env python
import logging
import argparse
import re


def append_fit_results():

    txt = args.fitlog.read()
    ts, te = re.findall("t=.*(\d+).*to.*(\d+)", txt)[0]
    mass, masserror = re.findall("mass.*(\d+\.\d+).*(\d+\.\d+)", txt)[0]
    amp, amperror = re.findall("amp.*(\d+\.\d+).*(\d+\.\d+)", txt)[0]
    qual = re.findall("Qual (\d+\.\d+)", txt)[0]
    fitcomment = "# fit({},{}) m={} e={} qual:{}".format(ts, te, mass, masserror, qual)
    print fitcomment,
    print args.filetoappend.read(),

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="append fit result (to an emass file)")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="increase output verbosity")
    parser.add_argument('fitlog', type=argparse.FileType('r'))
    parser.add_argument('filetoappend', type=argparse.FileType('r'))
    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.DEBUG)
        logging.debug("Verbose debuging mode activated")
    else:
        logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.INFO)

    append_fit_results()
