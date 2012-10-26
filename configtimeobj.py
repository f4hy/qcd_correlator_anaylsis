import numpy as np
import math
import logging


class Cfgtimeobj(object):
    """ A class to handle all the objects that are indexed as
    x[config][time] and provide some methods to simplify accessing the
    data in nonstandard ways, such as a time slice of all configs.
    """
    data = None
    sums = None
    average = None

    def __init__(self, datadict):
        self.data = datadict

        self.configs = datadict.keys()
        self.times = datadict[self.configs[0]].keys()
        self.numconfigs = len(self.configs)
        self.numtimes = len(self.times)

        logging.debug("created cto with configs %d, times %d", self.numconfigs, self.numtimes)

        dataitem = self[self.configs[0]][self.times[0]]

        self.datatype = type(dataitem)
        self.scalar = np.isscalar(dataitem)

        if __debug__:
            self.verify()

    @classmethod
    def fromDataDict(cls, datadict):
        return cls(datadict)

    @classmethod
    def fromListTuple(cls, listtuple):
        configs = list(range(len(listtuple)))

        data = {}
        for cfg in configs:
            rawcfgdata = listtuple[cfg]
            d = {}
            for rawtimedata in rawcfgdata:
                t = rawtimedata[0]  # First element is the time
                arraydata = np.asarray(tuple(rawtimedata)[1:])
                d[t] = arraydata  # rest is data e.g. real,imag
            data[cfg] = d
        return cls(data)

    def verify(self):

        if not __debug__:
            return True

        if not self.data:
            raise ValueError("data obejct empty or false")

        sizes = [len(v) for v in self.data.values()]  # map(len, self.data.values())
        if (sizes.count(sizes[0]) != len(sizes)):
            raise ValueError("Object size is inconsistant")

        for cfg in self.configs:
            for time in self.times:
                if type(self.data[cfg][time]) != self.datatype:
                    logging.error(str(type(self.data[cfg][time])))
                    raise TypeError("Not all data is the same type")
                if self.data[cfg][time] is None:
                    raise ValueError("indexed value is none")

        logging.debug("Cfg Time Object verified for consistancy")
        return True

    def __getitem__(self, key):
        return self.data[key]

    def __setitem__(self, key, value):
        self.data[key] = value

    def get(self, config=None, time=None):

        if config is None and time is None:
            return self.data
        if time is None:
            return self.data[config]
        if config is None:
            return {cfg: v[time] for cfg, v in self.data.iteritems()}
        else:
            return self.data[config][time]

    def compatible(self, otherobj):

        if self.configs == otherobj.configs and self.times == otherobj.times:
            logging.debug("verified objects are compatible")
            return True
        else:
            return False

    def indexes(self):
        return (self.configs, self.times)

    def average_over_times(self):
        return {config: math.fsum(self.get(config=config).values()) / float(self.numtimes)
                for config in self.configs}

    def average_all(self):
        if not self.average:
            aot = self.average_over_times()
            self.average = math.fsum(aot.values()) / float(self.numconfigs)
        return self.average

    def average_over_configs(self):
        sums = self.sum_over_configs()
        return {time: sums[time] / float(self.numconfigs) for time in self.times}

    def sum_over_configs(self):
        if not self.sums:
            self.sums = {t: math.fsum(self.get(time=t).values())
                         for t in self.times}
        return self.sums

    def jackknifed_averages(self):
        sums = self.sum_over_configs()
        return {cfg: {t: (sums[t] - self.get(time=t, config=cfg)) / (self.numconfigs - 1)
                      for t in self.times} for cfg in self.configs}

    def jackknifed_full_average(self):
        total = self.average_all()
        N = float(self.numconfigs)
        Njk = float(self.numconfigs - 1)
        return {cfg: (N * total - single) / Njk
                for cfg, single in self.average_over_times().iteritems()}

    def writefullfile(self, filename):
        outfile = open(filename, 'w')
        for cfg in self.configs:
            for t in self.times:
                if self.datatype is type(np.array(1.0)):
                    outfile.write("{!r},   {!r}\n".format(t, self.get(config=cfg, time=t)[0]))
                else:
                    outfile.write("{!r},   {!r}\n".format(t, self.get(config=cfg, time=t)))
        outfile.close()

    def writeeachconfig(self, filename):
        for config in self.configs:
            outfile = open(filename + '.' + str(config), 'w')
            for time in self.times:
                outfile.write("{!r},   {!r}\n".format(time, self.get(config=config, time=time)))
            outfile.close()
