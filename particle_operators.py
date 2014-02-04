#!/usr/bin/env python
import argparse
import ConfigParser
import logging
import readinput
import sys

momentum_map = {0: "AR", 1: "OA", 2: "PD", 3: "CD", 4: "OA", 5: "PSQ5"}
default_file = sys.path[0]+"/particles.ini"


class particleDatabase():

    def __init__(self, datafile=default_file):
        self.datafile = datafile
        logging.debug("initializing particleDatabase object")
        self.config = ConfigParser.SafeConfigParser()
        logging.debug("Reading database from {}".format(self.datafile))
        self.config.read(self.datafile)

    def read_op(self, particle, irrep, momentum):
        logging.debug("reading {}_{} from database".format(particle, momentum))
        key = "{}_{}".format(particle, momentum_map[momentum])
        try:
            return self.config.get(key, irrep)
        except ConfigParser.NoOptionError:
            logging.warn("No operator found for {} {}".format(key, irrep))
            sucsess = self.add_op_entry(key, irrep)
            if sucsess:
                logging.warn("Should have updated database, Re-trying")
                return self.read_op(particle, irrep, momentum)
            else:
                return None
        except ConfigParser.NoSectionError:
            logging.warn("No particle found for {}".format(key))
            self.config.add_section(key)
            sucsess = self.add_op_entry(key, irrep)
            if sucsess:
                logging.warn("Should have updated database, Re-trying")
                return self.read_op(particle, irrep, momentum)
            else:
                return None

    def add_op_entry(self, key, irrep):
        logging.info("please enter operator choice for {} {}".format(key, irrep))
        op = readinput.askoperator("{} {}".format(key, irrep))
        print "# user input operator", op
        if op:
            self.config.set(key, irrep, op)
            with open(self.datafile, "wb") as configfile:
                self.config.write(configfile)
            return True
        return False


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Simple interface for accessing particle database entries")
    parser.add_argument('particle', metavar='particle_name', type=str, help='particle name')
    parser.add_argument('momentum', metavar='momentum', type=int, help='particle momentum')
    parser.add_argument('irrep', metavar='irrep', type=str, help='particle irrep')
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="increase output verbosity")

    args = parser.parse_args()
    if args.verbose:
        logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.DEBUG)
        logging.debug("Verbose debuging mode activated")
    else:
        logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.INFO)
    mydatabase = particleDatabase()
    print mydatabase.read_op(args.particle, args.irrep, args.momentum)
