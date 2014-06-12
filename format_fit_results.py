#!/usr/bin/env python
import logging
import argparse
import re


def format_fit_results(filewild):

    level = 0
    while True:
        try:
            filename = filewild.format(level)
            logging.debug("opening file {}".format(filename))
            with open(filewild.format(level)) as infile:

                txt = infile.read()
                #print txt
                mass, masserror = re.findall("mass .*?(\d+\.\d+).*?(\d+\.\d+)", txt)[0]
                amp, amperror = re.findall("amp .*?(\d+\.\d+).*?(\d+\.\d+)", txt)[0]
                print "{}, {}, {}, {}, {}".format(level, amp, amperror, mass, masserror)

        except IOError:
            logging.debug("No file {}".format(filename))
            break
        level += 1

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="combine fit results to many levels into a single file")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="increase output verbosity")
    parser.add_argument('filewild', metavar='f', type=str, help='wildcard for fit outputs')
    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.DEBUG)
        logging.debug("Verbose debuging mode activated")
    else:
        logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.INFO)

    format_fit_results(args.filewild)
