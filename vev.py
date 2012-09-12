import math
import logging


class Vev(object):
    """Vev object
    """

    def __init__(self, datadict):
        """ Creates a vev object, needs a dict with value at each config"""
        self.data = datadict
        self.configs = datadict.keys()
        self.numconfigs = len(datadict.keys())
        logging.debug("created vev with %d configs", self.numconfigs)

    def __getitem__(self, key):
        return self.data[key]

    def __setitem__(self, key, value):
        self.data[key] = value

    def get(self, config):
        return self.data[config]

    def average(self):
        return math.fsum(self.data.values()) / float(self.numconfigs)

    def jackknife(self):
        return {cfg: (math.fsum(self.data.values()) - self.get(cfg)) / (self.numconfigs - 1)
                for cfg in self.configs}

    def writefullfile(self, filename):
        outfile = open(filename, 'w')
        for config in self.configs:
            outfile.write("%f\n" % (self[config]))
        outfile.close()
