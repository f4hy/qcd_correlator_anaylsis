#!/usr/bin/env python2
import argparse
import logging
import itertools


psqr_map = {0: [0,0,0], 1: [0,0,1], 2: [0,1,1], 3: [1,1,1], 4: [0,0,2],
           5: [0,1,2], 6: [1,1,2], 8: [0,2,2], 9: [0,0,3], -9: [1,2,2]}


def all_permutations(psqr, outputformat="{}{}{}", nobackies=False):
    if psqr == 9:
        logging.warn("There are two ways of making psqr=9, 003 and 122, use -9 for the 122 case")

    reference = psqr_map[psqr]
    all = set(itertools.permutations(reference))
    for i in itertools.product([1,-1], repeat=3):
        flipped = [a*b for a,b in zip(i,reference)]
        all.update(set(itertools.permutations(flipped)))

    if nobackies:
        for i in sorted(list(all), key=sum, reverse=True):
            if i not in all:
                continue
            reverse = tuple(-1*e for e in i)
            print i, reverse
            all.discard(reverse)

    joined = [outputformat.format(*mom) for mom in all]
    return joined



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="average data files")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="increase output verbosity")
    parser.add_argument("-o", "--output_stub", type=str, required=False,
                        help="stub of name to write output to")
    parser.add_argument("-p", "--psqr", type=int, required=True,
                    help="input momentum squared, produce all permuations")
    parser.add_argument("-nb", "--nobackies", action="store_true", required=False,
                        help="don't put in reverses'")
    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.DEBUG)
        logging.debug("Verbose debuging mode activated")
    else:
        logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.INFO)

    print " ".join(all_permutations(args.psqr, nobackies=args.nobackies))
