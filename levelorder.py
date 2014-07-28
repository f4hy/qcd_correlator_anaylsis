#!/usr/bin/env python
import logging
import argparse
import re


def print_results(results):
    for r in results:
        level, amp, amperror, mass, masserror, comment = r
        print "{}, {}, {}, {}, {}, {}".format(level, amp, amperror, mass, masserror, comment)

def print_order(results):
    s = [x[0] for x in sorted(results, key=lambda t: float(t[3]))]
    for r in s:
        print r

def format_fit_results(filewild):

    results = []
    level = 0
    while True:
        try:
            filename = filewild.format(level)
            logging.debug("opening file {}".format(filename))
            with open(filewild.format(level)) as infile:
                comment = ""
                txt = infile.read()
                if txt == "":
                    logging.debug("empty log file!")
                    level+=1
                    continue
                if "mass2" not in txt:
                    comment = " # single exp fit"
                mass, masserror = re.findall("mass .*?(\d+\.\d+).*?(\d+\.\d+)", txt)[0]
                amp, amperror = re.findall("amp .*?(\d+\.\d+).*?(\d+\.\d+)", txt)[0]
                results.append((level, amp, amperror, mass, masserror, comment))

        except IOError:
            logging.debug("No file {}".format(filename))
            break
        level += 1
    #print_results(results)
    print_order(results)


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
