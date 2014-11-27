#!/usr/bin/env python2
import sys
import re
import logging
import argparse

def format_op(op):
    "format operator from ruby format to david format"
    try:
        base, rep, p1,c1,o1, p2,c2,o2 = re.match('@oplist.push\("(.*)"\)', op).group(1).split()
    except ValueError:
        return

    iso, n1, n2 =  base.split("_")
    p1 = "".join(p1.replace("[P=(", "").replace(")", "").split(","))
    p2 = "".join(p2.replace("[P=(", "").replace(")", "").split(","))
    o1=o1.strip(']\"').replace("_","")
    o2=o2.strip(']\"').replace("_","")
    rep=rep.replace("_","")

    print "{}{}-{}p{}{}{}-{}p{}{}{}".format(iso, rep, n1, p1, o1, c1, n2, p2, o2, c2)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="translate ruby format to colin format")
    logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.DEBUG)
    parser.add_argument('strings', metavar='s', type=str, nargs='+',
                        help='expected level to determine opname from')
    args = parser.parse_args()

    if args.strings == ["-"]:
        args.strings = []
        for line in sys.stdin:
            args.strings.append( line.strip())

    for op in args.strings:
        format_op(op)
